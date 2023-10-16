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
//! \file MultiphaseBoundaryConditionsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseBoundaryConditionsBlockVisitor.h"
#include "BC.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "BC.h"
#include "Block3D.h"
#include "BCArray3D.h"
#include "LBMKernel.h"

MultiphaseBoundaryConditionsBlockVisitor::MultiphaseBoundaryConditionsBlockVisitor() :
Block3DVisitor(0, D3Q27System::MAXLEVEL)
{

}
//////////////////////////////////////////////////////////////////////////
MultiphaseBoundaryConditionsBlockVisitor::~MultiphaseBoundaryConditionsBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void MultiphaseBoundaryConditionsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if (block->getRank() == grid->getRank())
   {
      SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());

      if (!kernel)
      {
         throw UbException(UB_EXARGS, "LBMKernel in " + block->toString() + "is not exist!");
      }

      SPtr<BCSet> bcSet = kernel->getBCSet();

      if (!bcSet)
      {
         throw UbException(UB_EXARGS,"Boundary Conditions Processor is not exist!" );
      }

      SPtr<BCArray3D> bcArray = bcSet->getBCArray();

      bool compressible = kernel->getCompressible();
      real collFactorL = kernel->getCollisionFactorL();
	  real collFactorG = kernel->getCollisionFactorG();
	  real collFactorPh = 1.0/kernel->getPhaseFieldRelaxation();
	  real densityRatio = kernel->getDensityRatio();
	  real phiL = kernel->getPhiL();
	  real phiH = kernel->getPhiH();
      //int level = block->getLevel();

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;
      int maxX1 = (int)bcArray->getNX1();
      int maxX2 = (int)bcArray->getNX2();
      int maxX3 = (int)bcArray->getNX3();
      SPtr<BoundaryConditions> bcPtr;

      bcSet->clearBC();

      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
	  SPtr<DistributionArray3D> distributionsH = kernel->getDataSet()->getHdistributions();
      SPtr<DistributionArray3D> distributionsH2 = kernel->getDataSet()->getH2distributions();

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
                     char alg = bcPtr->getBCStrategyType();
                     SPtr<BCStrategy> bca = bcMap[alg];
                     
                     if (bca)
                     {
                        bca = bca->clone();
                        bca->setNodeIndex(x1, x2, x3);
                        bca->setBcPointer(bcPtr);
                        //bca->addDistributions(distributions, distributionsH);
						bca->addDistributions(distributions);
						bca->addDistributionsH(distributionsH);
                        if (distributionsH2)
                            bca->addDistributionsH2(distributionsH2);
                        bca->setCompressible(compressible);
                        bca->setBcArray(bcArray);
                        bcSet->addBC(bca);
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseBoundaryConditionsBlockVisitor::addBC(SPtr<BC> bc)
{
   bcMap.insert(std::make_pair(bc->getBCStrategyType(), bc->getBCStrategy()));
}



