#include "BoundaryConditionsBlockVisitor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "Grid3DSystem.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "ThinWallNoSlipBCAlgorithm.h"

BoundaryConditionsBlockVisitor::BoundaryConditionsBlockVisitor() :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
{

}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionsBlockVisitor::~BoundaryConditionsBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionsBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if (block->getRank() == grid->getRank())
   {
      LBMKernelPtr kernel = block->getKernel();
      BCProcessorPtr bcProcessor = kernel->getBCProcessor();

      if (!bcProcessor)
      {
         throw UbException(UB_EXARGS,"Boundary Conditions Processor is not exist!" );
      }

      BCArray3D& bcArray = bcProcessor->getBCArray();

      bool compressible = kernel->getCompressible();
      double collFactor = kernel->getCollisionFactor();
      int level = block->getLevel();

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;
      int maxX1 = (int)bcArray.getNX1();
      int maxX2 = (int)bcArray.getNX2();
      int maxX3 = (int)bcArray.getNX3();
      BoundaryConditionsPtr bcPtr;

      bcProcessor->clearBC();

      DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();

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
                     char alg = bcPtr->getBcAlgorithmType();
                     BCAlgorithmPtr bca = bcMap[alg];
                     
                     if (bca)
                     {
                        bca = bca->clone();
                        bca->addNode(x1, x2, x3);
                        bca->addBcPointer(bcPtr);
                        bca->addDistributions(distributions);
                        bca->setCollFactor(collFactor);
                        bca->setCompressible(compressible);
                        bcProcessor->addBC(bca);
                        bca->setBcArray(BCArray3DPtr(&bcArray));
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionsBlockVisitor::addBC(BCAdapterPtr bc)
{
   bcMap.insert(std::make_pair(bc->getBcAlgorithmType(), bc->getAlgorithm()));
}



