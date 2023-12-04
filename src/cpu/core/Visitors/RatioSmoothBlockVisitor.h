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
//! \file RatioSmoothBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef RatioSmoothBlockVisitor_H
#define RatioSmoothBlockVisitor_H

#include <string>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

class RatioSmoothBlockVisitor : public Block3DVisitor
{
public:
    RatioSmoothBlockVisitor(int levelDepth, bool includeNotActiveBlocks = false);

    ~RatioSmoothBlockVisitor() override = default;

    bool expandsByAdaptation() { return this->expandBlocks; }

    void setExpandByAdaptation(bool expandBlocks);

    int getLevelRatio() { return this->maxLevelRatio; }
    bool isIterative() { return true; }

    void setLevelRatio(int ratio);

    int getStartLevel();
    int getStopLevel();

    void setStartLevel(int level);
    void setStopLevel(int level);

    std::string getSpecificDescription();

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

protected:
    bool lookForExpand(SPtr<Grid3D> grid, const int &ix1, const int &ix2, const int &ix3, const int &level);
    bool lookForCollapse(SPtr<Grid3D> grid, const int &ix1, const int &ix2, const int &ix3, const int &level);

private:
    int maxLevelRatio;
    bool expandBlocks;
    int levelDepth;
    bool includeNotActiveBlocks;
};

#endif // RatioSmoothBlockVisitor_H
