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
//! \addtogroup gpu_grid_tests grid
//! \ingroup gpu_GridGenerator_tests GridGenerator
//! \{
//=======================================================================================
# include <gmock/gmock.h>
#include <string>
# include "Distribution.h"

#include "grid/distributions/D3Q27.h"
#include "lbm/constants/D3Q27.h"
using namespace vf::lbm::dir;

TEST(DistributionTest, DistributionReturnsCorrectDirections)
{
    Distribution dist = DistributionHelper::getDistribution27();

    EXPECT_THAT(dist.directions[dP00][0], testing::Eq(DIR_27_E_X));
    EXPECT_THAT(dist.directions[dP00][1], testing::Eq(DIR_27_E_Y));
    EXPECT_THAT(dist.directions[dP00][2], testing::Eq(DIR_27_E_Z));
    EXPECT_THAT(dist.dirs[dP00 * 3    ], testing::Eq(DIR_27_E_X));
    EXPECT_THAT(dist.dirs[dP00 * 3 + 1], testing::Eq(DIR_27_E_Y));
    EXPECT_THAT(dist.dirs[dP00 * 3 + 2], testing::Eq(DIR_27_E_Z));

    EXPECT_THAT(dist.directions[d00M][0], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.directions[d00M][1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.directions[d00M][2], testing::Eq(DIR_27_B_Z));
    EXPECT_THAT(dist.dirs[d00M * 3    ], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.dirs[d00M * 3 + 1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.dirs[d00M * 3 + 2], testing::Eq(DIR_27_B_Z));
    
    EXPECT_THAT(dist.directions[d000][0], testing::Eq(0));
    EXPECT_THAT(dist.directions[d000][1], testing::Eq(0));
    EXPECT_THAT(dist.directions[d000][2], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3    ], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3 + 1], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3 + 2], testing::Eq(0));

    EXPECT_THAT(dist.directions[dPP0][0], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.directions[dPP0][1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.directions[dPP0][2], testing::Eq(DIR_27_NE_Z));
    EXPECT_THAT(dist.dirs[dPP0 * 3    ], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.dirs[dPP0 * 3 + 1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.dirs[dPP0 * 3 + 2], testing::Eq(DIR_27_NE_Z));

    EXPECT_THAT(dist.directions[d0MP][0], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.directions[d0MP][1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.directions[d0MP][2], testing::Eq(DIR_27_TS_Z));
    EXPECT_THAT(dist.dirs[d0MP * 3    ], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.dirs[d0MP * 3 + 1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.dirs[d0MP * 3 + 2], testing::Eq(DIR_27_TS_Z));

    EXPECT_THAT(dist.directions[dPPP][0], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.directions[dPPP][1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.directions[dPPP][2], testing::Eq(DIR_27_TNE_Z));
    EXPECT_THAT(dist.dirs[dPPP * 3    ], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.dirs[dPPP * 3 + 1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.dirs[dPPP * 3 + 2], testing::Eq(DIR_27_TNE_Z));

    EXPECT_THAT(dist.directions[dMMM][0], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.directions[dMMM][1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.directions[dMMM][2], testing::Eq(DIR_27_BSW_Z));
    EXPECT_THAT(dist.dirs[dMMM * 3    ], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.dirs[dMMM * 3 + 1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.dirs[dMMM * 3 + 2], testing::Eq(DIR_27_BSW_Z));
}
//! \}
