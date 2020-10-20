#if defined VF_METIS && defined VF_MPI

#include "PePartitioningGridVisitor.h"
#include "Block3D.h"
#include "Communicator.h"
#include "CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "UbLogger.h"
#include <math.h>
#include <shared_mutex>

#include "DemCoProcessor.h"

using namespace std;

PePartitioningGridVisitor::PePartitioningGridVisitor(SPtr<Communicator> comm, std::shared_ptr<DemCoProcessor> dem)
    : Grid3DVisitor(), comm(comm), dem(dem)
{
    forest = dynamicPointerCast<PePhysicsEngineSolverAdapter>(dem->getPhysicsEngineSolver())->getForest();
}
//////////////////////////////////////////////////////////////////////////
PePartitioningGridVisitor::~PePartitioningGridVisitor() {}
//////////////////////////////////////////////////////////////////////////
void PePartitioningGridVisitor::visit(SPtr<Grid3D> grid)
{
    UBLOG(logDEBUG1, "PePartitioningGridVisitor::visit() - start");

    collectData(grid);
    distributePartitionData(grid);

    UBLOG(logDEBUG1, "PePartitioningGridVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void PePartitioningGridVisitor::collectData(SPtr<Grid3D> grid)
{
    // int minInitLevel = grid->getCoarsestInitializedLevel();
    // int maxInitLevel = grid->getFinestInitializedLevel();

    walberla::uint_t peRank;

    for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
        forest->getProcessRank(peRank, blockIt->getId());
        vector<SPtr<Block3D>> blocks;
        walberla::AABB aabb = blockIt->getAABB();

        // getBlocksByCuboid((double)aabb.xMin(), (double)aabb.yMin(), (double)aabb.zMin(), (double)aabb.xMax(),
        // (double)aabb.yMax(), (double)aabb.zMax(), blocks, grid); for (SPtr<Block3D> block : blocks)
        //{
        //   ids.push_back(block->getGlobalID());
        //   ranks.push_back((int)peRank);
        //}
        SPtr<Block3D> block = getBlockByMinUniform((double)aabb.xMin(), (double)aabb.yMin(), (double)aabb.zMin(), grid);
        if (block) {
            ids.push_back(block->getGlobalID());
            ranks.push_back((int)peRank);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
// void PePartitioningGridVisitor::getBlocksByCuboid(double minX1, double minX2, double minX3, double maxX1, double
// maxX2, double maxX3, std::vector<SPtr<Block3D>>& blocks, SPtr<Grid3D> grid)
//{
//   int coarsestLevel = grid->getCoarsestInitializedLevel();
//   int finestLevel   = grid->getFinestInitializedLevel();
//
//   SPtr<CoordinateTransformation3D> trafo = grid->getCoordinateTransformator();
//
//   //////////////////////////////////////////////////////////////////////////
//   //MINIMALE BLOCK-INDIZES BESTIMMEN
//   //
//   //min:
//   double dMinX1 = trafo->transformForwardToX1Coordinate(minX1, minX2, minX3)*(1<<finestLevel);
//   double dMinX2 = trafo->transformForwardToX2Coordinate(minX1, minX2, minX3)*(1<<finestLevel);
//   double dMinX3 = trafo->transformForwardToX3Coordinate(minX1, minX2, minX3)*(1<<finestLevel);
//
//   //Achtung, wenn minX1 genau auf grenze zwischen zwei bloecken -> der "kleinere" muss genommen werden,
//   //da beim Transformieren der "groessere" Index rauskommt
//   int iMinX1 = (int)dMinX1; //if (UbMath::zero(dMinX1-iMinX1)) iMinX1-=1;
//   int iMinX2 = (int)dMinX2; //if (UbMath::zero(dMinX2-iMinX2)) iMinX2-=1;
//   int iMinX3 = (int)dMinX3; //if (UbMath::zero(dMinX3-iMinX3)) iMinX3-=1;
//
//   //max (hier kann die Zusatzabfrage vernachlaessigt werden):
//   int iMaxX1 = (int)(trafo->transformForwardToX1Coordinate(maxX1, maxX2, maxX3)*(1<<finestLevel));
//   int iMaxX2 = (int)(trafo->transformForwardToX2Coordinate(maxX1, maxX2, maxX3)*(1<<finestLevel));
//   int iMaxX3 = (int)(trafo->transformForwardToX3Coordinate(maxX1, maxX2, maxX3)*(1<<finestLevel));
//
//   SPtr<Block3D> block;
//
//   //set, um doppelte bloecke zu vermeiden, die u.U. bei periodic auftreten koennen
//   std::set<SPtr<Block3D>> blockset;
//   for (int level=coarsestLevel; level<=finestLevel; level++)
//   {
//      //damit bei negativen werten auch der "kleinere" genommen wird -> floor!
//      int minx1 = (int)std::floor((double)iMinX1/(1<<(finestLevel-level)));
//      int minx2 = (int)std::floor((double)iMinX2/(1<<(finestLevel-level)));
//      int minx3 = (int)std::floor((double)iMinX3/(1<<(finestLevel-level)));
//
//      int maxx1 = iMaxX1/(1<<(finestLevel-level));
//      int maxx2 = iMaxX2/(1<<(finestLevel-level));
//      int maxx3 = iMaxX3/(1<<(finestLevel-level));
//
//      for (int ix1=minx1; ix1<maxx1; ix1++)
//         for (int ix2=minx2; ix2<maxx2; ix2++)
//            for (int ix3=minx3; ix3<maxx3; ix3++)
//               if ((block=grid->getBlock(ix1, ix2, ix3, level)))
//               {
//                  blockset.insert(block);
//               }
//   }
//
//   blocks.resize(blockset.size());
//   std::copy(blockset.begin(), blockset.end(), blocks.begin());
//}

SPtr<Block3D> PePartitioningGridVisitor::getBlockByMinUniform(double minX1, double minX2, double minX3,
                                                              SPtr<Grid3D> grid)
{
    SPtr<CoordinateTransformation3D> trafo = grid->getCoordinateTransformator();

    int ix1 = (int)trafo->transformForwardToX1Coordinate(minX1, minX2, minX3);
    int ix2 = (int)trafo->transformForwardToX2Coordinate(minX1, minX2, minX3);
    int ix3 = (int)trafo->transformForwardToX3Coordinate(minX1, minX2, minX3);

    return grid->getBlock(ix1, ix2, ix3, 0);
}

//////////////////////////////////////////////////////////////////////////
void PePartitioningGridVisitor::distributePartitionData(SPtr<Grid3D> grid)
{
    std::vector<int> totalIDs;
    std::vector<int> totalRanks;

    assert(ids.size() != 0);
    assert(ranks.size() != 0);

    comm->allGather(ids, totalIDs);
    comm->allGather(ranks, totalRanks);

    assert(totalIDs.size() == totalRanks.size());
    for (int i = 0; i < totalIDs.size(); i++) {
        SPtr<Block3D> block = grid->getBlock(totalIDs[i]);
        if (block)
            block->setRank(totalRanks[i]);
    }
}
//////////////////////////////////////////////////////////////////////////

#endif // VF_METIS
