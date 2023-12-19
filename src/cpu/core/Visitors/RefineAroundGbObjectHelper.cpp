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
#include "RefineAroundGbObjectHelper.h"
#include <parallel/Communicator.h>
#include "OverlapBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include <D3Q27System.h>
#include <D3Q27TriFaceMeshInteractor.h>
#include <Grid3D.h>

RefineAroundGbObjectHelper::RefineAroundGbObjectHelper(SPtr<Grid3D> grid, int refineLevel,
                                                       SPtr<D3Q27TriFaceMeshInteractor> objectIter,
                                                       real startDistance, real stopDistance,
                                                       std::shared_ptr<vf::parallel::Communicator> comm)
    : grid(grid), refineLevel(refineLevel), objectIter(objectIter), startDistance(startDistance),
      stopDistance(stopDistance), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
RefineAroundGbObjectHelper::~RefineAroundGbObjectHelper(void) = default;
//////////////////////////////////////////////////////////////////////////
void RefineAroundGbObjectHelper::refine()
{
    using namespace vf::lbm::dir;

    UBLOG(logDEBUG5, "RefineCrossAndInsideGbObjectHelper: refine - start");

    int rank = grid->getRank();
    grid->setRank(0);

    objectIter->refineBlockGridToLevel(refineLevel, startDistance, stopDistance);

    RatioBlockVisitor ratioVisitor(refineLevel);
    grid->accept(ratioVisitor);

    RatioSmoothBlockVisitor ratioSmoothVisitor(refineLevel);
    grid->accept(ratioSmoothVisitor);

    OverlapBlockVisitor overlapVisitor(refineLevel, false);
    grid->accept(overlapVisitor);

    std::vector<int> dirs;
    for (int i = (int)dP00; i <= (int)d0MP; i++) {
        dirs.push_back(i);
    }
    SetInterpolationDirsBlockVisitor interDirsVisitor(dirs);
    grid->accept(interDirsVisitor);

    grid->setRank(rank);

    UBLOG(logDEBUG5, "RefineCrossAndInsideGbObjectHelper: refine - end");
}

//! \}
