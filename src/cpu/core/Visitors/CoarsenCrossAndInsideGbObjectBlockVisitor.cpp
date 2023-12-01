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
//! \file CoarsenCrossAndInsideGbObjectBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================
#include "CoarsenCrossAndInsideGbObjectBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include <geometry3d/GbObject3D.h>

CoarsenCrossAndInsideGbObjectBlockVisitor::CoarsenCrossAndInsideGbObjectBlockVisitor() : Block3DVisitor() {}
//////////////////////////////////////////////////////////////////////////
CoarsenCrossAndInsideGbObjectBlockVisitor::CoarsenCrossAndInsideGbObjectBlockVisitor(SPtr<GbObject3D> geoObject,
                                                                                     int fineLevel, int coarseLevel)
    : Block3DVisitor(fineLevel, fineLevel), geoObject(geoObject), coarseLevel(coarseLevel)
{
}
//////////////////////////////////////////////////////////////////////////
CoarsenCrossAndInsideGbObjectBlockVisitor::~CoarsenCrossAndInsideGbObjectBlockVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void CoarsenCrossAndInsideGbObjectBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    int fineLevel = block->getLevel();
    if (notActive && block->isNotActive())
        return;
    if (fineLevel > this->getStopLevel())
        return;

    UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
    UbTupleDouble3 deltas = grid->getBlockLengths(block);
    if (geoObject->isCellInsideOrCuttingGbObject3D(val<1>(coords), val<2>(coords), val<3>(coords),
                                                   val<1>(coords) + val<1>(deltas), val<2>(coords) + val<2>(deltas),
                                                   val<3>(coords) + val<3>(deltas))) {
        grid->collapseBlock(block->getX1(), block->getX2(), block->getX3(), fineLevel, coarseLevel);
    }

    return;
}
//////////////////////////////////////////////////////////////////////////
