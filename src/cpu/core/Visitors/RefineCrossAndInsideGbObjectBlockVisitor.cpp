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
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"

#include "Block3D.h"
#include "Grid3D.h"
#include <geometry3d/GbObject3D.h>

RefineCrossAndInsideGbObjectBlockVisitor::RefineCrossAndInsideGbObjectBlockVisitor() : Block3DVisitor() {}
//////////////////////////////////////////////////////////////////////////
RefineCrossAndInsideGbObjectBlockVisitor::RefineCrossAndInsideGbObjectBlockVisitor(SPtr<GbObject3D> geoObject,
                                                                                   int refineLevel)
    : Block3DVisitor(0, refineLevel - 1), geoObject(geoObject)
{
}
//////////////////////////////////////////////////////////////////////////
RefineCrossAndInsideGbObjectBlockVisitor::~RefineCrossAndInsideGbObjectBlockVisitor() = default;
//////////////////////////////////////////////////////////////////////////
void RefineCrossAndInsideGbObjectBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    int level = block->getLevel();
    if (notActive && block->isNotActive())
        return;
    if (level > this->getStopLevel())
        return;

    UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
    UbTupleDouble3 deltas = grid->getBlockLengths(block);
    if (geoObject->isCellInsideOrCuttingGbObject3D(val<1>(coords), val<2>(coords), val<3>(coords),
                                                   val<1>(coords) + val<1>(deltas), val<2>(coords) + val<2>(deltas),
                                                   val<3>(coords) + val<3>(deltas))) {
        grid->expandBlock(block->getX1(), block->getX2(), block->getX3(), level);
    }

    return;
}
//////////////////////////////////////////////////////////////////////////

//! \}
