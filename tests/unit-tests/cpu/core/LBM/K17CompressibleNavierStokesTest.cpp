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
//! \addtogroup cpu_Connectors_tests Connectors
//! \ingroup cpu_core_tests core
//! \{

#include <gtest/gtest.h>
//#include <gmock/gmock.h>
#include <array>
#include <memory>

#include <cpu/core/LBM/K17CompressibleNavierStokes.h>

class K17CompressibleNavierStokesTest : public ::testing::Test {
protected:
    K17CompressibleNavierStokes solver;

};

TEST_F(K17CompressibleNavierStokesTest, ValidLimiterIsSetCorrectly) {
    std::array<real, 3> input = {0.5, 0.5, 0.5};
    solver.setQuadricLimiter(input);
    EXPECT_EQ(solver.getQuadricLimiter(), input);
}

TEST_F(K17CompressibleNavierStokesTest, InvalidLimiterDefaultsToFallback) {
    std::array<real, 3> invalidInput = {-10.0, 0.0, 9999.0};
    solver.setQuadricLimiter(invalidInput);

    std::array<real, 3> expected = {0.01, 0.01, 0.01};

    EXPECT_EQ(solver.getQuadricLimiter(), expected);
}

TEST_F(K17CompressibleNavierStokesTest, InvalidInitOfLimiterDefaultsToFallback) {
    std::array<real, 3> invalidInit;
    solver.setQuadricLimiter(invalidInit);

    std::array<real, 3> expected = {0.01, 0.01, 0.01};

    EXPECT_EQ(solver.getQuadricLimiter(), expected);
}

//! \}