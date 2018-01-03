#include "SpongeLayerBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

using namespace std;

SpongeLayerBlockVisitor::SpongeLayerBlockVisitor(GbCuboid3DPtr boundingBox, LBMReal collFactor) : 
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL), 
   boundingBox(boundingBox),
   collFactor(collFactor)
{

}
//////////////////////////////////////////////////////////////////////////
SpongeLayerBlockVisitor::~SpongeLayerBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if (block->getRank() == grid->getRank())
   {
      UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
      UbTupleDouble3 blockLengths = grid->getBlockLengths(block);

      double minX1 = val<1>(org);
      double minX2 = val<2>(org);
      double minX3 = val<3>(org);
      double maxX1 = val<1>(org)+val<1>(blockLengths);
      double maxX2 = val<2>(org)+val<2>(blockLengths);
      double maxX3 = val<3>(org)+val<3>(blockLengths);

      if (boundingBox->isCellInsideGbObject3D(minX1, minX2, minX3, maxX1, maxX2, maxX3))
      {
         LBMKernelPtr kernel = block->getKernel();
         //kernel->setCollisionFactor(LBMSystem::calcCollisionFactor(0.01, block->getLevel()));
         kernel->setCollisionFactor(collFactor);
      }
   }
}


