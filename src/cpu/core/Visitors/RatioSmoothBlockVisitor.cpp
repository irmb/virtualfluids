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
//! \file RatioSmoothBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include "RatioSmoothBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"

RatioSmoothBlockVisitor::RatioSmoothBlockVisitor(int levelDepth, bool includeNotActiveBlocks)
    : Block3DVisitor(D3Q27System::MAXLEVEL, 0), maxLevelRatio(1), expandBlocks(true), levelDepth(levelDepth),
      includeNotActiveBlocks(includeNotActiveBlocks)
{
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
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
bool RatioSmoothBlockVisitor::lookForExpand(SPtr<Grid3D> grid, const int &ix1, const int &ix2, const int &ix3,
                                            const int &level)
{
    std::vector<SPtr<Block3D>> neighbors;
    grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
    int nix1, nix2, nix3, nlev;
    for (size_t i = 0; i < neighbors.size(); i++) {
        if ((neighbors[i]->isActive() || includeNotActiveBlocks) && neighbors[i]->getLevel() > level) {
            nix1 = (neighbors)[i]->getX1();
            nix2 = (neighbors)[i]->getX2();
            nix3 = (neighbors)[i]->getX3();
            nlev = (neighbors)[i]->getLevel();

            std::vector<SPtr<Block3D>> neighbors1;
            grid->getAllNeighbors(nix1, nix2, nix3, nlev, nlev + 1, neighbors1);
            for (size_t j = 0; j < neighbors1.size(); j++) {
                if ((neighbors1[j]->isActive() || includeNotActiveBlocks) &&
                    neighbors1[j]->getLevel() > level + this->maxLevelRatio) {
                    return true;
                }
            }
        }
    }
    return false;
}
//////////////////////////////////////////////////////////////////////////
bool RatioSmoothBlockVisitor::lookForCollapse(SPtr<Grid3D> grid, const int &ix1, const int &ix2, const int &ix3,
                                              const int &level)
{
    std::vector<SPtr<Block3D>> neighbors;
    grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
    for (size_t i = 0; i < neighbors.size(); i++) {
        if ((neighbors[i]->isActive() || includeNotActiveBlocks) &&
            neighbors[i]->getLevel() < level - this->maxLevelRatio) {
            throw UbException(UB_EXARGS, " not implemented till now");
            return true;
        }
    }

    return false;
}
//////////////////////////////////////////////////////////////////////////
void RatioSmoothBlockVisitor::setExpandByAdaptation(bool expandBlocks)
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
void RatioSmoothBlockVisitor::setLevelRatio(int ratio)
{
    if (ratio < 1)
        throw UbException(UB_EXARGS, "illegal ratio specified");
    this->maxLevelRatio = ratio;
}
//////////////////////////////////////////////////////////////////////////
std::string RatioSmoothBlockVisitor::getSpecificDescription()
{
    std::string str("Ratio:");
    return str;
}
//////////////////////////////////////////////////////////////////////////
int RatioSmoothBlockVisitor::getStartLevel()
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
int RatioSmoothBlockVisitor::getStopLevel()
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
void RatioSmoothBlockVisitor::setStartLevel(int level)
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
void RatioSmoothBlockVisitor::setStopLevel(int level)
{
    if (this->expandBlocks) {
        if (level <= Block3DVisitor::getStartLevel())
            Block3DVisitor::setStopLevel(level);
    } else {
        if (level >= Block3DVisitor::getStartLevel())
            Block3DVisitor::setStopLevel(level);
    }
}
