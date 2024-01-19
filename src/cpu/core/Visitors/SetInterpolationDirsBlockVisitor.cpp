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
#include "SetInterpolationDirsBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include <D3Q27System.h>

SetInterpolationDirsBlockVisitor::SetInterpolationDirsBlockVisitor(std::vector<int> &dirs)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), dirs(dirs)
{
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationDirsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

    int ix1, ix2, ix3, level;
    ix1   = block->getX1();
    ix2   = block->getX2();
    ix3   = block->getX3();
    level = block->getLevel();
    using namespace D3Q27System;
    if (level == 0)
        return;

    SPtr<Block3D> parentblock = grid->getSuperBlock(ix1, ix2, ix3, level);
    if (!parentblock)
        return;

    for (int dir : dirs) {
        SPtr<Block3D> nblock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);
        if (!nblock) {
            SPtr<Block3D> p_nblock = grid->getNeighborBlock(dir, parentblock);

            if (p_nblock) {
                bool flagDir;
                switch (dir) {
                    case dPP0:
                        checkFlagDir(grid, dP00, d0P0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dMM0:
                        checkFlagDir(grid, dM00, d0M0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dPM0:
                        checkFlagDir(grid, dP00, d0M0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dMP0:
                        checkFlagDir(grid, dM00, d0P0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dP0P:
                        checkFlagDir(grid, dP00, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dM0M:
                        checkFlagDir(grid, dM00, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dP0M:
                        checkFlagDir(grid, dP00, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dM0P:
                        checkFlagDir(grid, dM00, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case d0PP:
                        checkFlagDir(grid, d0P0, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case d0MM:
                        checkFlagDir(grid, d0M0, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case d0PM:
                        checkFlagDir(grid, d0P0, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case d0MP:
                        checkFlagDir(grid, d0M0, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dPPP:
                        checkFlagDir(grid, dP00, d0P0, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dMMP:
                        checkFlagDir(grid, dM00, d0M0, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dPMP:
                        checkFlagDir(grid, dP00, d0M0, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dMPP:
                        checkFlagDir(grid, dM00, d0P0, d00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dPPM:
                        checkFlagDir(grid, dP00, d0P0, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dMMM:
                        checkFlagDir(grid, dM00, d0M0, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dPMM:
                        checkFlagDir(grid, dP00, d0M0, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case dMPM:
                        checkFlagDir(grid, dM00, d0P0, d00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                }

                block->setInterpolationFlagFC(dir);
                parentblock->setInterpolationFlagCF(dir);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
// void SetInterpolationDirsBlockVisitor::checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, bool &flagDirection, int
// ix1, int ix2, int ix3, int level)
//{
//   SPtr<Block3D> block1 = grid->getNeighborBlock(dir1, ix1, ix2, ix3, level);
//   SPtr<Block3D> block2 = grid->getNeighborBlock(dir2, ix1, ix2, ix3, level);
//   if (!((block1 && block2)  ||  (!block1 && !block2)))
//      flagDirection = false;
//   else
//      flagDirection = true;
//}

void SetInterpolationDirsBlockVisitor::checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, bool &flagDirection, int ix1,
                                                    int ix2, int ix3, int level)
{
    SPtr<Block3D> block1 = grid->getNeighborBlock(dir1, ix1, ix2, ix3, level);
    SPtr<Block3D> block2 = grid->getNeighborBlock(dir2, ix1, ix2, ix3, level);

    SPtr<Block3D> pblock  = grid->getSuperBlock(ix1, ix2, ix3, level);
    SPtr<Block3D> pblock1 = grid->getNeighborBlock(dir1, pblock);
    SPtr<Block3D> pblock2 = grid->getNeighborBlock(dir2, pblock);

    if (!((block1 && block2) || (!block1 && !block2)) || !((pblock1 && pblock2) || (!pblock1 && !pblock2)))
        flagDirection = false;
    else
        flagDirection = true;
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationDirsBlockVisitor::checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, int dir3,
                                                    bool &flagDirection, int ix1, int ix2, int ix3, int level)
{
    SPtr<Block3D> block1 = grid->getNeighborBlock(dir1, ix1, ix2, ix3, level);
    SPtr<Block3D> block2 = grid->getNeighborBlock(dir2, ix1, ix2, ix3, level);
    SPtr<Block3D> block3 = grid->getNeighborBlock(dir3, ix1, ix2, ix3, level);
    if (!((block1 && block2 && block3) || (!block1 && !block2 && !block3)))
        flagDirection = false;
    else
        flagDirection = true;
}

//! \}
