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
#include "SetForcingBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "LBMSystem.h"

SetForcingBlockVisitor::SetForcingBlockVisitor(real forcingX1, real forcingX2, real forcingX3)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), forcingX1(forcingX1), forcingX2(forcingX2), forcingX3(forcingX3)
{
    ftype = 0;
}
//////////////////////////////////////////////////////////////////////////
SetForcingBlockVisitor::SetForcingBlockVisitor(const mu::Parser &muForcingX1, const mu::Parser &muForcingX2,
                                               const mu::Parser &muForcingX3)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), muForcingX1(muForcingX1), muForcingX2(muForcingX2),
      muForcingX3(muForcingX3)

{
    ftype = 1;
}
//////////////////////////////////////////////////////////////////////////
SetForcingBlockVisitor::SetForcingBlockVisitor(const std::string &sForcingX1, const std::string &sForcingX2,
                                               const std::string &sForcingX3)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), sForcingX1(sForcingX1), sForcingX2(sForcingX2), sForcingX3(sForcingX3)

{
    ftype = 2;
}
//////////////////////////////////////////////////////////////////////////
void SetForcingBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (block->getRank() == grid->getRank()) {
        SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
        if (!kernel)
            throw UbException(UB_EXARGS, "LBMKernel is not exist");

        switch (ftype) {
            case 0:
                kernel->setForcingX1(forcingX1);
                kernel->setForcingX2(forcingX2);
                kernel->setForcingX3(forcingX3);
                kernel->setWithForcing(true);
                break;
            case 1:
                kernel->setForcingX1(muForcingX1);
                kernel->setForcingX2(muForcingX2);
                kernel->setForcingX3(muForcingX3);
                kernel->setWithForcing(true);
                break;
            case 2:
                kernel->setForcingX1(sForcingX1);
                kernel->setForcingX2(sForcingX2);
                kernel->setForcingX3(sForcingX3);
                kernel->setWithForcing(true);
                break;
            default:
                kernel->setForcingX1(0.0);
                kernel->setForcingX2(0.0);
                kernel->setForcingX3(0.0);
                kernel->setWithForcing(false);
                break;
        }
    }
}

//! \}
