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
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include "RatioBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"

RatioBlockVisitor::RatioBlockVisitor(int levelDepth, bool includeNotActiveBlocks)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), maxLevelRatio(1), expandBlocks(true), levelDepth(levelDepth),
      includeNotActiveBlocks(includeNotActiveBlocks)
{
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    int ix1, ix2, ix3, level;
    ix1   = block->getX1();
    ix2   = block->getX2();
    ix3   = block->getX3();
    level = block->getLevel();

    if (block->isActive() || includeNotActiveBlocks) {
        if (this->expandBlocks) {
            if (this->lookForExpand(grid, ix1, ix2, ix3, level)) {
                grid->expandBlock(ix1, ix2, ix3, level);
            }
        } else {
            if (this->lookForCollapse(grid, ix1, ix2, ix3, level)) {
                grid->collapseBlock(ix1, ix2, ix3, level, levelDepth);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
bool RatioBlockVisitor::lookForExpand(SPtr<Grid3D> grid, const int &ix1, const int &ix2, const int &ix3,
                                      const int &level)
{
    std::vector<SPtr<Block3D>> neighbors;
    grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
    for (size_t i = 0; i < neighbors.size(); i++) {
        if ((neighbors[i]->isActive() || includeNotActiveBlocks) &&
            neighbors[i]->getLevel() > level + this->maxLevelRatio) {
            return true;
        }
    }

    return false;
}
//////////////////////////////////////////////////////////////////////////
bool RatioBlockVisitor::lookForCollapse(SPtr<Grid3D> grid, const int &ix1, const int &ix2, const int &ix3,
                                        const int &level)
{
    std::vector<SPtr<Block3D>> neighbors;
    grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
    for (size_t i = 0; i < neighbors.size(); i++) {
        if ((neighbors[i]->isActive() || includeNotActiveBlocks) &&
            neighbors[i]->getLevel() < level - this->maxLevelRatio) {
            return true;
        }
    }

    return false;
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setExpandByAdaptation(bool expandBlocks)
{
    if (this->expandBlocks != expandBlocks) {
        this->expandBlocks = expandBlocks;

        int l1 = Block3DVisitor::getStartLevel();
        int l2 = Block3DVisitor::getStopLevel();

        if (expandBlocks) {
            if (l1 < l2) {
                Block3DVisitor::setStartLevel(l2);
                Block3DVisitor::setStopLevel(l1);
            }
        } else {
            if (l2 < l1) {
                Block3DVisitor::setStartLevel(l2);
                Block3DVisitor::setStopLevel(l1);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setLevelRatio(int ratio)
{
    if (ratio < 1)
        throw UbException(UB_EXARGS, "illegal ratio specified");
    this->maxLevelRatio = ratio;
}
//////////////////////////////////////////////////////////////////////////
std::string RatioBlockVisitor::getSpecificDescription()
{
    std::string str("Ratio:");
    return str;
}
//////////////////////////////////////////////////////////////////////////
int RatioBlockVisitor::getStartLevel()
{
    int l1 = Block3DVisitor::getStartLevel();
    int l2 = Block3DVisitor::getStopLevel();

    if (this->expandBlocks) {
        if (l2 < l1)
            return (l1);
        else
            return (l2);
    } else {
        if (l2 < l1)
            return (l2);
        else
            return (l1);
    }
}
//////////////////////////////////////////////////////////////////////////
int RatioBlockVisitor::getStopLevel()
{
    int l1 = Block3DVisitor::getStartLevel();
    int l2 = Block3DVisitor::getStopLevel();

    if (this->expandBlocks) {
        if (l2 < l1)
            return (l2);
        else
            return (l1);
    } else {
        if (l2 < l1)
            return (l1);
        else
            return (l2);
    }
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setStartLevel(int level)
{
    if (this->expandBlocks) {
        if (level >= Block3DVisitor::getStopLevel())
            Block3DVisitor::setStartLevel(level);
    } else {
        if (level <= Block3DVisitor::getStopLevel())
            Block3DVisitor::setStartLevel(level);
    }
}
//////////////////////////////////////////////////////////////////////////
void RatioBlockVisitor::setStopLevel(int level)
{
    if (this->expandBlocks) {
        if (level <= Block3DVisitor::getStartLevel())
            Block3DVisitor::setStopLevel(level);
    } else {
        if (level >= Block3DVisitor::getStartLevel())
            Block3DVisitor::setStopLevel(level);
    }
}

//! \}
