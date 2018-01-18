#include "SpongeLayerBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"

#include "Grid3D.h"
#include "Block3D.h"
#include "ILBMKernel.h"
#include "UbException.h"

#include "CompressibleCumulant4thOrderViscosityLBMKernel.h"

#include <numerics/geometry3d/GbCuboid3D.h>

using namespace std;

SpongeLayerBlockVisitor::SpongeLayerBlockVisitor() : Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
{
   
}
//////////////////////////////////////////////////////////////////////////
SpongeLayerBlockVisitor::~SpongeLayerBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if (!boundingBox)
   {
      UB_THROW(UbException(UB_EXARGS, "The bounding box isn't set!"));
   }
   if (!kernel)
   {
      UB_THROW(UbException(UB_EXARGS, "The kernel isn't set!"));
   }
   if (kernel && (block->getRank() == grid->getRank()))
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
         LBMReal collFactor = block->getKernel()->getCollisionFactor();
         kernel->setCollisionFactor(collFactor);
         kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
         kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
         kernel->setBlock(block);
         SPtr<LBMKernel> newKernel = kernel->clone();

         SPtr<DataSet3D> dataSet = block->getKernel()->getDataSet();
         if (!dataSet)
         {
            UB_THROW(UbException(UB_EXARGS, "It is not possible to change a DataSet in kernel! Old DataSet is not exist!"));
         }

         newKernel->setDataSet(dataSet);

         SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
         if (!bcProc)
         {
            UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
         }
         newKernel->setBCProcessor(bcProc);

         double oldCollFactor = newKernel->getCollisionFactor();

         int ibX1 = block->getX1();

         UbTupleInt3 ixMin = grid->getBlockIndexes(boundingBox->getX1Minimum(),boundingBox->getX2Minimum(),boundingBox->getX3Minimum());
         UbTupleInt3 ixMax = grid->getBlockIndexes(boundingBox->getX1Maximum(),boundingBox->getX2Maximum(),boundingBox->getX3Maximum());

         int ibMax = val<1>(ixMax)-val<1>(ixMin)+1;
         double index = (double)(ibX1-val<1>(ixMin)+1);

         double newCollFactor = oldCollFactor - (oldCollFactor-1.0)/(double)(ibMax)*index;

         newKernel->setCollisionFactor(newCollFactor);
         block->setKernel(newKernel);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::setBoundingBox(SPtr<GbCuboid3D> bb)
{
   boundingBox = bb;
}
//////////////////////////////////////////////////////////////////////////
void SpongeLayerBlockVisitor::setKernel(SPtr<LBMKernel> k)
{
   kernel = k;
}


