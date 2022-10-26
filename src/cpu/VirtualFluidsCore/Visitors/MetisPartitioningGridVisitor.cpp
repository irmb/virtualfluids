#if defined VF_METIS && defined VF_MPI

#include "MetisPartitioningGridVisitor.h"
#include "Block3D.h"
#include <mpi/Communicator.h>
#include "D3Q27System.h"
#include "Grid3D.h"
#include <cmath>

using namespace std;

MetisPartitioningGridVisitor::MetisPartitioningGridVisitor(std::shared_ptr<vf::mpi::Communicator> comm, GraphType graphType, int numOfDirs,
                                                           MetisPartitioner::PartType partType, bool threads,
                                                           int numberOfThreads)
    : Grid3DVisitor(), numberOfThreads(numberOfThreads), numOfDirs(numOfDirs), comm(comm), threads(threads),
      graphType(graphType), partType(partType)
{
    numberOfProcesses = comm->getNumberOfProcesses();
}
//////////////////////////////////////////////////////////////////////////
MetisPartitioningGridVisitor::~MetisPartitioningGridVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::visit(SPtr<Grid3D> grid)
{
    UBLOG(logDEBUG1, "MetisPartitioningGridVisitor::visit() - start");

    this->clear();

    bundleRoot      = comm->getBundleRoot();
    bundleID        = comm->getBundleID();
    numberOfBundles = comm->getNumberOfBundles();
    if (numberOfBundles > 1) {
        if (bundleRoot == bundleID && processRoot == processID)
            collectData(grid, numberOfBundles, BUNDLE);
        comm->broadcast(blockID);
        comm->broadcast(parts);
        distributePartitionData(grid, BUNDLE);
        this->clear();
    }

    processRoot = comm->getProcessRoot();
    processID   = comm->getProcessID();
    /*int numberOfProcesses = comm->getNumberOfProcesses();*/
    if (numberOfProcesses > 1) {
        int temp = bundleID;
        for (int i = 0; i < numberOfBundles; i++) {
            if (bundleRoot == bundleID && processRoot == processID) {
                bundleID = i;
                // numberOfProcesses = comm->getNumberOfProcessesInBundle(i);
                collectData(grid, numberOfProcesses, PROCESS);
                bundleID = temp;
            }
            comm->broadcast(blockID);
            // UBLOG(logINFO, "blockID="<<blockID.size());
            comm->broadcast(parts);
            // UBLOG(logINFO, "parts="<<parts.size());
            distributePartitionData(grid, PROCESS);
        }
    }

    if (threads) {
        if (numberOfThreads > 1) {
            collectData(grid, numberOfThreads, THREAD);
            distributePartitionData(grid, THREAD);
        }
    }
    UBLOG(logDEBUG1, "MetisPartitioningGridVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::collectData(SPtr<Grid3D> grid, int nofSegments, PartLevel level)
{
    clear();

    switch (graphType) {
        case LevelIntersected:
            buildMetisGraphLevelIntersected(grid, nofSegments, level);
            break;
        case LevelBased:
            buildMetisGraphLevelBased(grid, nofSegments, level);
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::distributePartitionData(SPtr<Grid3D> grid, PartLevel level)
{
    SPtr<Block3D> block;

    for (size_t p = 0; p < parts.size(); p++) {
        block = grid->getBlock(blockID[p]);
        if (block) {
            switch (level) {
                case BUNDLE:
                    block->setBundle(parts[p]);
                    break;
                case PROCESS:
                    if (numberOfBundles == 1) {
                        block->setRank(parts[p]);
                    } else {
                        block->setLocalRank(parts[p]);
                        block->setRank(comm->getProcessID(block->getBundle(), parts[p]));
                    }
                    break;
                case THREAD:
                    block->setPart(parts[p]);
                    break;
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::buildMetisGraphLevelIntersected(SPtr<Grid3D> grid, int nofSegments, PartLevel level)
{
    int edges                       = 0;
    const int edgeWeight            = 1;
    const int edgeWeightChildFactor = 8;
    int n                           = 0;

    for (Grid3D::BlockIDMap::value_type b : grid->getBlockIDs()) {
        SPtr<Block3D> block = b.second;
        if (this->getPartitionCondition(block, level)) {
            block->setLocalID(n);
            blockID.push_back(block->getGlobalID());
            n++;
        }
    }

    MetisPartitioner metis;

    for (Grid3D::BlockIDMap::value_type b : grid->getBlockIDs()) {
        const SPtr<Block3D> block = b.second;
        if (this->getPartitionCondition(block, level)) {
            metis.xadj.push_back(edges);
            // the weights of the vertices are 2^level of grid (1, 2, 4, 8 .....) 1<<level
            metis.vwgt.push_back((idx_t)(1 << block->getLevel()));

            for (int dir = 0; dir <= numOfDirs; dir++) {
                SPtr<Block3D> neighBlock = grid->getNeighborBlock(dir, block);
                if (neighBlock) {
                    if (this->getPartitionCondition(neighBlock, level)) {
                        edges++;
                        metis.adjwgt.push_back(edgeWeight);
                        metis.adjncy.push_back(neighBlock->getLocalID());
                    }
                }
            }
            vector<SPtr<Block3D>> subBlocks;
            grid->getSubBlocks(block, 1, subBlocks);
            for (SPtr<Block3D> subBlock : subBlocks) {
                if (subBlock) {
                    if (this->getPartitionCondition(subBlock, level)) {
                        edges++;
                        metis.adjwgt.push_back(edgeWeight * edgeWeightChildFactor);
                        metis.adjncy.push_back(subBlock->getLocalID());
                    }
                }
            }
        }
    }

    metis.xadj.push_back(static_cast<idx_t>(metis.adjncy.size()));
    if ((metis.adjncy.size() % 2) != 0)
        throw UbException(
            UB_EXARGS,
            "number of edges is odd - probable adjncy-vector doesn't contain all pairs (A->B) and (B->A)!!!");

    metis.partition(nofSegments, partType);
    parts = metis.part;
}
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::buildMetisGraphLevelBased(SPtr<Grid3D> grid, int nofSegments, PartLevel level)
{
    int minInitLevel = grid->getCoarsestInitializedLevel();
    int maxInitLevel = grid->getFinestInitializedLevel();

    for (int l = minInitLevel; l <= maxInitLevel; l++) {
        int n = 0;
        vector<SPtr<Block3D>> blockVector;
        grid->getBlocks(l, blockVector);
        vector<SPtr<Block3D>> tBlockID;

        for (SPtr<Block3D> block : blockVector) {
            if (this->getPartitionCondition(block, level)) {
                block->setLocalID(n);
                blockID.push_back(block->getGlobalID());
                tBlockID.push_back(block);
                n++;
            }
        }

        if (tBlockID.size() == 0) {
            UB_THROW(UbException(UB_EXARGS, "Blocks for decomposition don't exist!"));
        }

        MetisPartitioner metis;

        const int vertexWeight = 1;
        int edges              = 0;

        for (SPtr<Block3D> block : tBlockID) {
            metis.xadj.push_back(edges);
            metis.vwgt.push_back(vertexWeight);

            for (int dir = 0; dir <= numOfDirs; dir++) {
                SPtr<Block3D> neighBlock = grid->getNeighborBlock(dir, block);
                if (neighBlock) {
                    if (this->getPartitionCondition(neighBlock, level)) {
                        edges++;
                        metis.adjwgt.push_back(getEdgeWeight(dir));
                        metis.adjncy.push_back(neighBlock->getLocalID());
                    }
                }
            }
        }
        metis.xadj.push_back(static_cast<idx_t>(metis.adjncy.size()));
        if ((metis.adjncy.size() % 2) != 0)
            throw UbException(
                UB_EXARGS,
                "number of edges is odd - probable adjncy-vector doesn't contain all pairs (A->B) and (B->A)!!!");

        int nofBlocks    = grid->getNumberOfBlocks(l);
        int tnofSegments = nofSegments;
        if (nofBlocks < nofSegments) {
            tnofSegments = nofBlocks;
        }
        metis.partition(tnofSegments, partType);

        for (idx_t p : metis.part) {
            parts.push_back(p);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
bool MetisPartitioningGridVisitor::getPartitionCondition(SPtr<Block3D> block, PartLevel level)
{
    if (level == BUNDLE) {
        return true;
    } else if (level == PROCESS) {
        if (block->getBundle() == bundleID) {
            return true;
        }
    } else if (level == THREAD) {
        if (block->getBundle() == bundleID && block->getRank() == processID) {
            return true;
        }
    }

    return false;
}
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::clear()
{
    blockID.clear();
    parts.clear();
}
//////////////////////////////////////////////////////////////////////////
int MetisPartitioningGridVisitor::getEdgeWeight(int dir)
{
    using namespace D3Q27System;
    if (dir <= DIR_00M) {
        return 100;
    } else if (dir >= DIR_PP0 && dir <= DIR_0MP) {
        return 10;
    } else if (dir >= DIR_PPP) {
        return 1;
    }

    //    FIXME: non-void function does not return a value in all control paths
    return 0;
}
//////////////////////////////////////////////////////////////////////////
void MetisPartitioningGridVisitor::setNumberOfProcesses(int np) { numberOfProcesses = np; }

#endif // VF_METIS
