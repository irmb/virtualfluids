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
#include "RefineCrossAndInsideGbObjectHelper.h"
#include "CheckRatioBlockVisitor.h"
#include <parallel/Communicator.h>
#include "OverlapBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include <D3Q27System.h>
#include <GbObject3D.h>
#include <Grid3D.h>

RefineCrossAndInsideGbObjectHelper::RefineCrossAndInsideGbObjectHelper(SPtr<Grid3D> grid, int maxRefineLevel,
                                                                       std::shared_ptr<vf::parallel::Communicator> comm)
    : grid(grid), maxRefineLevel(maxRefineLevel), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
RefineCrossAndInsideGbObjectHelper::~RefineCrossAndInsideGbObjectHelper(void) = default;
//////////////////////////////////////////////////////////////////////////
void RefineCrossAndInsideGbObjectHelper::refine()
{
    UBLOG(logDEBUG5, "RefineCrossAndInsideGbObjectHelper: refine - start");

    if (comm->isRoot()) {
        int size = (int)objects.size();

        for (int i = 0; i < size; i++) {
            RefineCrossAndInsideGbObjectBlockVisitor refVisitor(objects[i], levels[i]);
            grid->accept(refVisitor);
        }

        // RatioBlockVisitor ratioVisitor(maxRefineLevel);
        // grid->accept(ratioVisitor);

        // RatioSmoothBlockVisitor ratioSmoothVisitor(maxRefineLevel);
        // grid->accept(ratioSmoothVisitor);

        RatioBlockVisitor ratioVisitor(maxRefineLevel);
        CheckRatioBlockVisitor checkRatio(maxRefineLevel);
        int count = 0;

        do {
            grid->accept(ratioVisitor);
            checkRatio.resetState();
            grid->accept(checkRatio);
            UBLOG(logINFO, "count = " << count++ << " state = " << checkRatio.getState());
        } while (!checkRatio.getState());

        OverlapBlockVisitor overlapVisitor(maxRefineLevel, false);
        grid->accept(overlapVisitor);
    }

    grid->updateDistributedBlocks(comm);

    std::vector<int> dirs;

    for (int i = D3Q27System::FSTARTDIR; i <= D3Q27System::FENDDIR; i++) {
        dirs.push_back(i);
    }
    SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
    grid->accept(interDirsVisitor);
    UBLOG(logDEBUG5, "RefineCrossAndInsideGbObjectHelper: refine - end");
}
//////////////////////////////////////////////////////////////////////////
void RefineCrossAndInsideGbObjectHelper::addGbObject(SPtr<GbObject3D> object, int refineLevel)
{
    objects.push_back(object);
    levels.push_back(refineLevel);
}

//! \}
