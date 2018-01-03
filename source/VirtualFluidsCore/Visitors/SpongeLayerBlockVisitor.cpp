#include "SpongeLayerBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

#include "Grid3D.h"
#include "Block3D.h"
#include "ILBMKernel.h"

#include <numerics/geometry3d/GbCuboid3D.h>

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
         ILBMKernelPtr kernel = block->getKernel();
         double oldCollFactor = kernel->getCollisionFactor();

         int ibX1 = block->getX1();

         UbTupleInt3 ixMin = grid->getBlockIndexes(boundingBox->getX1Minimum(),boundingBox->getX2Minimum(),boundingBox->getX3Minimum());
         UbTupleInt3 ixMax = grid->getBlockIndexes(boundingBox->getX1Maximum(),boundingBox->getX2Maximum(),boundingBox->getX3Maximum());

         int ibMax = val<1>(ixMax)-val<1>(ixMin)+1;
         double index = ibX1-val<1>(ixMin)+1;

         double newCollFactor = oldCollFactor - (oldCollFactor-1.0)/ibMax*index;

         kernel->setCollisionFactor(newCollFactor);
      }
   }
}


