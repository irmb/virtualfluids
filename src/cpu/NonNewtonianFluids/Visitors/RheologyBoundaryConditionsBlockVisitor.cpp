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
#include "ILBMKernel.h"
#include "BCStrategyType.h"

#include "NonNewtonianFluids/BoundaryConditions/ThixotropyDensityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyVelocityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyNonReflectingOutflowBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyVelocityWithDensityBCStrategy.h"


BoundaryConditionsBlockVisitor::BoundaryConditionsBlockVisitor() : Block3DVisitor(0, D3Q27System::MAXLEVEL)
{
}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionsBlockVisitor::~BoundaryConditionsBlockVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void BoundaryConditionsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (block->getRank() == grid->getRank()) {
        SPtr<ILBMKernel> kernel = block->getKernel();

        if (!kernel) {
            throw UbException(UB_EXARGS, "LBMKernel in " + block->toString() + "is not exist!");
        }

        SPtr<BCSet> bcSet = kernel->getBCSet();

        if (!bcSet) {
            throw UbException(UB_EXARGS, "Boundary Conditions Processor is not exist!");
        }

        SPtr<BCArray3D> bcArray = bcSet->getBCArray();

        bool compressible = kernel->getCompressible();
        real collFactor = kernel->getCollisionFactor();

        int minX1 = 0;
        int minX2 = 0;
        int minX3 = 0;
        int maxX1 = (int)bcArray->getNX1();
        int maxX2 = (int)bcArray->getNX2();
        int maxX3 = (int)bcArray->getNX3();
        SPtr<BoundaryConditions> bcPtr;

        bcSet->clearBC();

        SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

        for (int x3 = minX3; x3 < maxX3; x3++) {
            for (int x2 = minX2; x2 < maxX2; x2++) {
                for (int x1 = minX1; x1 < maxX1; x1++) {
                    if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
                        if ((bcPtr = bcArray->getBC(x1, x2, x3)) != NULL) {
                            char alg              = bcPtr->getBCStrategyType();
                            SPtr<BCStrategy> bca = bcMap[alg];

                            if (bca) {
                                bca = bca->clone();
                                bca->setBlock(block);
                                bca->setNodeIndex(x1, x2, x3);
                                bca->setBcPointer(bcPtr);
                                bca->addDistributions(distributions);

                                if (alg == BCStrategyType::ThixotropyVelocityBCStrategy)
                                    std::static_pointer_cast<ThixotropyVelocityBCStrategy>(bca)->addDistributionsH(
                                        kernel->getDataSet()->getHdistributions());
                                if (alg == BCStrategyType::ThixotropyDensityBCStrategy)
                                    std::static_pointer_cast<ThixotropyDensityBCStrategy>(bca)->addDistributionsH(
                                        kernel->getDataSet()->getHdistributions());
                                if (alg == BCStrategyType::ThixotropyNoSlipBCStrategy)
                                    std::static_pointer_cast<ThixotropyNoSlipBCStrategy>(bca)->addDistributionsH(
                                        kernel->getDataSet()->getHdistributions());
                                if (alg == BCStrategyType::ThixotropyNonReflectingOutflowBCStrategy)
                                    std::static_pointer_cast<ThixotropyNonReflectingOutflowBCStrategy>(bca)
                                        ->addDistributionsH(kernel->getDataSet()->getHdistributions());
                                if (alg == BCStrategyType::ThixotropyVelocityWithDensityBCStrategy)
                                    std::static_pointer_cast<ThixotropyVelocityWithDensityBCStrategy>(bca)
                                        ->addDistributionsH(kernel->getDataSet()->getHdistributions());

                                bca->setCollFactor(collFactor);
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
void BoundaryConditionsBlockVisitor::addBC(SPtr<BC> bc)
{
    bcMap.insert(std::make_pair(bc->getBCStrategyType(), bc->getBCStrategy()));
}
