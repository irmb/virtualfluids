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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "BoundaryConditionsBlockVisitor.h"
#include "BC.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "EsoSplit.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "BC.h"
#include "Block3D.h"
#include "BCArray3D.h"
#include "ILBMKernel.h"

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
                            char bcStrategyKey = bcPtr->getBCStrategyKey();
                            SPtr<BCStrategy> bca = BCStrategyRegistry::getInstance()->getBCStrategy(bcStrategyKey);

                            if (bca) {
                                bca = bca->clone();
                                bca->setBlock(block);
                                bca->setNodeIndex(x1, x2, x3);
                                bca->setBcPointer(bcPtr);
                                bca->addDistributions(distributions);
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


//! \}
