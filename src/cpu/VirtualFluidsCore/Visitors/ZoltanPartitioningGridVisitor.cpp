#if defined VF_ZOLTAN && defined VF_MPI

#include "ZoltanPartitioningGridVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "UbException.h"
#include "UbLogger.h"
#include <vector>

using namespace std;

ZoltanPartitioningGridVisitor::ZoltanPartitioningGridVisitor(std::shared_ptr<vf::mpi::Communicator> comm, int numOfDirs,
                                                             int numOfLocalParts)
    : comm(comm), numOfDirs(numOfDirs), numOfLocalParts(numOfLocalParts)
{
}
//////////////////////////////////////////////////////////////////////////
ZoltanPartitioningGridVisitor::~ZoltanPartitioningGridVisitor() {}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioningGridVisitor::visit(SPtr<Grid3D> grid)
{
    UBLOG(logDEBUG5, "ZoltanPartitioningPatchVisitor::visit() - start");

    // MPI_Comm mpi_comm = *((MPI_Comm*) comm->getNativeCommunicator());

    // ZoltanPartitioner zp(mpi_comm, comm->getProcessID(), numOfLocalParts);
    //
    // graph = zp.getGraphData();

    collectData(grid);

    // zp.partition();

    // repartGrid(grid, zp);

    UBLOG(logDEBUG5, "ZoltanPartitioningPatchVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioningGridVisitor::collectData(SPtr<Grid3D> grid)
{
    int myRank    = comm->getProcessID();
    int numOfProc = comm->getNumberOfProcesses();

    if (numOfProc < 2) {
        return;
    }

    int minInitLevel = grid->getCoarsestInitializedLevel();
    int maxInitLevel = grid->getFinestInitializedLevel();

    for (int l = minInitLevel; l <= maxInitLevel; l++) {
        MPI_Comm mpi_comm = *((MPI_Comm *)comm->getNativeCommunicator());
        ZoltanPartitioner zp(mpi_comm, comm->getProcessID(), numOfLocalParts);
        graph = zp.getGraphData();

        int n = 0;
        vector<SPtr<Block3D>> blockVector;
        grid->getBlocks(l, blockVector);

        if (blockVector.size() == 0) {
            UB_THROW(UbException(UB_EXARGS, "Blocks for decomposition don't exist!"));
        }

        // Verteilung von Ranks
        int rank = 0;
        for (SPtr<Block3D> block : blockVector) {
            block->setRank(rank);
            block->setPart(rank);
            rank++;
            if (rank > numOfProc - 1)
                rank = 0;
        }

        int vertices = 0;

        for (SPtr<Block3D> block : blockVector) {
            if (block->getRank() == myRank) {
                vertices++;

                vertexGID.push_back(block->getGlobalID());

                int edges = 0;
                for (int dir = 0; dir <= numOfDirs; dir++) {
                    SPtr<Block3D> neighBlock =
                        (grid->getNeighborBlock(dir, block->getX1(), block->getX2(), block->getX3(), l));

                    if (neighBlock) {
                        edges++;
                        nborGID.push_back(neighBlock->getGlobalID());
                        nborProc.push_back(neighBlock->getRank());
                    }
                }
                numEdges.push_back(edges);
            }
        }
        graph->numLocalVertices = vertices;
        graph->vnumEdges        = numEdges;
        graph->vvertexGID       = vertexGID;
        graph->vnborGID         = nborGID;
        graph->vnborProc        = nborProc;

        zp.partition();
        repartGrid(grid, zp);
    }
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioningGridVisitor::repartGrid(SPtr<Grid3D> grid, ZoltanPartitioner &zp)
{
    if (zp.areChanges()) {
        UBLOG(logDEBUG5, "ZoltanPartitioningPatchVisitor::repartGrid - start");
        vector<int> sExportGlobalGids, sExportToPart, sExportProcs;
        vector<int> rExportGlobalGids, rExportToPart, rExportProcs;

        zp.getExportData(sExportGlobalGids, sExportToPart, sExportProcs);

        comm->allGather(sExportGlobalGids, rExportGlobalGids);
        comm->allGather(sExportToPart, rExportToPart);
        comm->allGather(sExportProcs, rExportProcs);

        for (int i = 0; i < (int)rExportGlobalGids.size(); i++) {
            if (rExportGlobalGids[i] != -1) {
                SPtr<Block3D> block = grid->getBlock(rExportGlobalGids[i]);
                if (block) {
                    block->setRank(rExportProcs[i]);
                    block->setPart(rExportToPart[i]);
                }
            }
        }
        UBLOG(logDEBUG5, "ZoltanPartitioningPatchVisitor::repartGrid - end");
    }
}
#endif
