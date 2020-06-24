#include "BoundaryConditionsBlockVisitor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "Grid3DSystem.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "ThinWallNoSlipBCAlgorithm.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "BCAdapter.h"
#include "Block3D.h"
#include "BCArray3D.h"

BoundaryConditionsBlockVisitor::BoundaryConditionsBlockVisitor() :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
{

}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionsBlockVisitor::~BoundaryConditionsBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if (block->getRank() == grid->getRank())
   {
      SPtr<ILBMKernel> kernel = block->getKernel();

      if (!kernel)
      {
         throw UbException(UB_EXARGS, "LBMKernel in " + block->toString() + "is not exist!");
      }

      SPtr<BCProcessor> bcProcessor = kernel->getBCProcessor();

      if (!bcProcessor)
      {
         throw UbException(UB_EXARGS,"Boundary Conditions Processor is not exist!" );
      }

      SPtr<BCArray3D> bcArray = bcProcessor->getBCArray();

      bool compressible = kernel->getCompressible();
      double collFactor = kernel->getCollisionFactor();
      int level = block->getLevel();

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;
      int maxX1 = (int)bcArray->getNX1();
      int maxX2 = (int)bcArray->getNX2();
      int maxX3 = (int)bcArray->getNX3();
      SPtr<BoundaryConditions> bcPtr;

      bcProcessor->clearBC();

      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

      for (int x3 = minX3; x3 < maxX3; x3++)
      {
         for (int x2 = minX2; x2 < maxX2; x2++)
         {
            for (int x1 = minX1; x1 < maxX1; x1++)
            {
               if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3))
               {
                  if ((bcPtr = bcArray->getBC(x1, x2, x3)) != NULL)
                  {
                     char alg = bcPtr->getBcAlgorithmType();
                     SPtr<BCAlgorithm> bca = bcMap[alg];
                     
                     if (bca)
                     {
                        bca = bca->clone();
                        bca->setNodeIndex(x1, x2, x3);
                        bca->setBcPointer(bcPtr);
                        bca->addDistributions(distributions);
                        bca->setCollFactor(collFactor);
                        bca->setCompressible(compressible);
                        bca->setBcArray(bcArray);
                        bcProcessor->addBC(bca);
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionsBlockVisitor::addBC(SPtr<BCAdapter> bc)
{
   bcMap.insert(std::make_pair(bc->getBcAlgorithmType(), bc->getAlgorithm()));
}



