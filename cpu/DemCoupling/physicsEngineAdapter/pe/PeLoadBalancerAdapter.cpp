#include "PeLoadBalancerAdapter.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "CoordinateTransformation3D.h"
#include "UbLogger.h"

#include "core/debug/CheckFunctions.h"


PeLoadBalancerAdapter::PeLoadBalancerAdapter(SPtr<Grid3D> grid, unsigned numberOfProcesses, int rank) : grid(grid), numberOfProcesses(numberOfProcesses), rank(rank)
{

}

walberla::uint_t PeLoadBalancerAdapter::operator()(walberla::SetupBlockForest & forest, const walberla::uint_t numberOfProcesses, const walberla::memory_t perProcessMemoryLimit)
{
   std::vector< walberla::SetupBlock * > peBlocks;
   forest.getBlocks(peBlocks);

   for (auto peBlock = peBlocks.begin(); peBlock != peBlocks.end(); ++peBlock)
   {
      walberla::AABB aabb = (*peBlock)->getAABB();
      SPtr<Block3D> block = getBlockByMinUniform(aabb.xMin()+0.5*(aabb.xMax()-aabb.xMin()), aabb.yMin()+0.5*(aabb.yMax()-aabb.yMin()), aabb.zMin()+0.5*(aabb.zMax()-aabb.zMin()), grid);
      if (block)
      {
         (*peBlock)->assignTargetProcess((walberla::uint_t)block->getRank());
      }
      else
      {
         //TODO: the rank of pe blocks is not consistent with VF blocks 
         (*peBlock)->assignTargetProcess(0);
         //UBLOG(logINFO, "PeLoadBalancerAdapter::operator() peBlockId="<<(*peBlock)->getId());
      }
   }

   return numberOfProcesses;
}

SPtr<Block3D> PeLoadBalancerAdapter::getBlockByMinUniform(double minX1, double minX2, double minX3, SPtr<Grid3D> grid)
{
   SPtr<CoordinateTransformation3D> trafo = grid->getCoordinateTransformator();

   int ix1 = (int)trafo->transformForwardToX1Coordinate(minX1, minX2, minX3);
   int ix2 = (int)trafo->transformForwardToX2Coordinate(minX1, minX2, minX3);
   int ix3 = (int)trafo->transformForwardToX3Coordinate(minX1, minX2, minX3);

   return grid->getBlock(ix1, ix2, ix3, 0);
}
