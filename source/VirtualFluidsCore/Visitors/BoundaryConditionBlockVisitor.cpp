#include "BoundaryConditionBlockVisitor.h"
#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include "Grid3DSystem.h"

BoundaryConditionBlockVisitor::BoundaryConditionBlockVisitor() :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
{

}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionBlockVisitor::~BoundaryConditionBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if (block->getRank() == grid->getRank())
   {
      LBMKernel3DPtr kernel = block->getKernel();
      D3Q27ETBCProcessorPtr bcProcessor = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor());
      velocity = bcProcessor->getBC(BoundaryCondition::Velocity);
      density = bcProcessor->getBC(BoundaryCondition::Density);
      noSlip = bcProcessor->getBC(BoundaryCondition::NoSlip);
      slip = bcProcessor->getBC(BoundaryCondition::Slip);
      BCArray3D<D3Q27BoundaryCondition>& bcArray = bcProcessor->getBCArray();
      bool compressible = kernel->getCompressible();
      double collFactor = kernel->getCollisionFactor();
      int level = block->getLevel();

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;
      int maxX1 = (int)bcArray.getNX1();
      int maxX2 = (int)bcArray.getNX2();
      int maxX3 = (int)bcArray.getNX3();
      D3Q27BoundaryConditionPtr bcPtr;

      int vCount = 0;
      int dCount = 0;
      int nsCount = 0;
      int sCount = 0;

      if (velocity) velocity->clearData();
      if (density) density->clearData();
      if (noSlip)  noSlip->clearData();
      if (slip)    slip->clearData();

      EsoTwist3DPtr distributions = boost::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());

      for (int x3 = minX3; x3 < maxX3; x3++)
      {
         for (int x2 = minX2; x2 < maxX2; x2++)
         {
            for (int x1 = minX1; x1 < maxX1; x1++)
            {
               if (!bcArray.isSolid(x1, x2, x3) && !bcArray.isUndefined(x1, x2, x3))
               {
                  if ((bcPtr = bcArray.getBC(x1, x2, x3)) != NULL)
                  {
                     //velocity boundary condition
                     if (bcPtr->hasVelocityBoundary() && velocity)
                     {
                        velocity->addNode(x1, x2, x3);
                        velocity->addBcPointer(bcPtr);
                        vCount++;
                     }
                     //density boundary condition
                     else if (bcPtr->hasDensityBoundary() && density)
                     {
                        density->addNode(x1, x2, x3);
                        density->addBcPointer(bcPtr);
                        dCount++;
                     }
                     //no slip boundary condition
                     else if (bcPtr->hasNoSlipBoundary() && noSlip)
                     {
                        noSlip->addNode(x1, x2, x3);
                        noSlip->addBcPointer(bcPtr);
                        nsCount++;
                     }
                     //slip boundary condition
                     else if (bcPtr->hasSlipBoundary() && slip)
                     {
                        slip->addNode(x1, x2, x3);
                        slip->addBcPointer(bcPtr);
                        sCount++;
                     }
                  }
               }
            }
         }
      }
      if (vCount > 0)
      {
         velocity->addDistributions(distributions);
         velocity->setCollFactor(collFactor);
         velocity->setCompressible(compressible);
      }
      if (dCount > 0)
      {
         density->addDistributions(distributions);
         density->setCollFactor(collFactor);
         density->setCompressible(compressible);
      }
      if (nsCount > 0)
      {
         noSlip->addDistributions(distributions);
         noSlip->setCollFactor(collFactor);
         noSlip->setCompressible(compressible);
      }
      if (sCount > 0)
      {
         slip->addDistributions(distributions);
         slip->setCollFactor(collFactor);
         slip->setCompressible(compressible);
      }
   }
}


