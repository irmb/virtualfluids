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
#include "RefineInterGbObjectsVisitor.h"

#include "Block3D.h"
#include "Grid3D.h"
#include <geometry3d/GbObject3D.h>
#include "lbm/constants/D3Q27.h"

RefineInterGbObjectsBlockVisitor::RefineInterGbObjectsBlockVisitor() : Block3DVisitor(-1, -1) {}
//////////////////////////////////////////////////////////////////////////
RefineInterGbObjectsBlockVisitor::RefineInterGbObjectsBlockVisitor(SPtr<GbObject3D> includeGbObject3D,
                                                                   SPtr<GbObject3D> excludeGbObject3D, int startlevel,
                                                                   int stoplevel)
    : Block3DVisitor(startlevel, stoplevel)
{
    this->includeGbObjects3D.push_back(includeGbObject3D);
    this->excludeGbObjects3D.push_back(excludeGbObject3D);
}
//////////////////////////////////////////////////////////////////////////
RefineInterGbObjectsBlockVisitor::RefineInterGbObjectsBlockVisitor(std::vector<SPtr<GbObject3D>> includeGbObjects3D,
                                                                   std::vector<SPtr<GbObject3D>> excludeGbObjects3D,
                                                                   int startlevel, int stoplevel)
    : Block3DVisitor(startlevel, stoplevel)
{
    this->includeGbObjects3D = includeGbObjects3D;
    this->excludeGbObjects3D = excludeGbObjects3D;
}
//////////////////////////////////////////////////////////////////////////
void RefineInterGbObjectsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    UbTupleDouble3 coords = grid->getBlockWorldCoordinates(block);
    UbTupleDouble3 delta  = grid->getBlockLengths(block);

    real cellMinX1 = val<1>(coords);
    real cellMinX2 = val<2>(coords);
    real cellMinX3 = val<3>(coords);
    real cellMaxX1 = val<1>(coords) + val<1>(delta);
    real cellMaxX2 = val<2>(coords) + val<2>(delta);
    real cellMaxX3 = val<3>(coords) + val<3>(delta);

    bool insideInclude = false;
    for (size_t i = 0; i < includeGbObjects3D.size(); i++) {
        if (includeGbObjects3D[i]->isCellInsideOrCuttingGbObject3D(cellMinX1, cellMinX2, cellMinX3, cellMaxX1,
                                                                   cellMaxX2, cellMaxX3)) {
            insideInclude = true;
            break;
        }
    }

    bool insideExclude = false;
    for (size_t e = 0; e < excludeGbObjects3D.size(); e++) {
        if (excludeGbObjects3D[e]->isCellInsideGbObject3D(cellMinX1, cellMinX2, cellMinX3, cellMaxX1, cellMaxX2,
                                                          cellMaxX3)) {
            insideExclude = true;
            break;
        }
    }

    if (insideInclude && !insideExclude) {
        int ix1, ix2, ix3, level;
        ix1   = block->getX1();
        ix2   = block->getX2();
        ix3   = block->getX3();
        level = block->getLevel();
        grid->expandBlock(ix1, ix2, ix3, level);
    }
}

//! \}
