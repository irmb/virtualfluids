//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file BoundaryConditionsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#include "BoundaryConditionsBlockVisitor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "Grid3DSystem.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
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



