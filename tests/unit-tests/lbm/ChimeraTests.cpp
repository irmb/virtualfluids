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
//! \addtogroup lbm
//! \{
//! \author Soeren Peters
//=======================================================================================
#include <gmock/gmock.h>

#include <lbm/ChimeraTransformation.h>

#ifdef VF_DOUBLE_ACCURACY
#define REAL_EQ(a) testing::DoubleEq(a)
#else
#define REAL_EQ(a) testing::FloatEq(a)
#endif

/*
* InverseChimeraWithK
*/
TEST(ChimeraTest, forwardInverseChimeraWithK)
{
    real mfa = 1;
    real mfb = 1;
    real mfc = 1;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;
    const real Kinverse = 1 / K;

    vf::lbm::forwardInverseChimeraWithK(mfa, mfb, mfc, vv, v2, K, Kinverse);

    EXPECT_THAT(mfa, REAL_EQ(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, REAL_EQ(-4.)); // -(mfa + mfb + mfc + 1)
    EXPECT_THAT(mfc, REAL_EQ(6.));  // (mfa + mfc) + (mfa + mfb + mfc + 1)
}


TEST(ChimeraTest, backwardInverseChimeraWithK)
{
    // starting with the result values from the test above.
    real mfa = 3.;
    real mfb = -4.;
    real mfc = 6.;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;
    const real Kinverse = 1 / K;

    vf::lbm::backwardInverseChimeraWithK(mfa, mfb, mfc, vv, v2, K, Kinverse);

    // resulting in the start values from the test above.
    EXPECT_THAT(mfa, REAL_EQ(1.));
    EXPECT_THAT(mfb, REAL_EQ(1.));
    EXPECT_THAT(mfc, REAL_EQ(1.));
}

/*
* Chimera
*/
TEST(ChimeraTest, forwardChimera)
{
    real mfa = 1;
    real mfb = 1;
    real mfc = 1;

    const real vv = 1.;
    const real v2 = 1.;

    vf::lbm::forwardChimera(mfa, mfb, mfc, vv, v2);

    EXPECT_THAT(mfa, REAL_EQ(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, REAL_EQ(-3.)); // -(mfa + mfb + mfc)
    EXPECT_THAT(mfc, REAL_EQ(5.));  // (mfa + mfc) + (mfa + mfb + mfc)
}


TEST(ChimeraTest, backwardChimera)
{
    // starting with the result values from the test above.
    real mfa = 3.;
    real mfb = -3.;
    real mfc = 5.;

    const real vv = 1.;
    const real v2 = 1.;

    vf::lbm::backwardChimera(mfa, mfb, mfc, vv, v2);

    // resulting in the start values from the test above.
    EXPECT_THAT(mfa, REAL_EQ(1.));
    EXPECT_THAT(mfb, REAL_EQ(1.));
    EXPECT_THAT(mfc, REAL_EQ(1.));
}

/*
* ChimeraWithK
*/
TEST(ChimeraTest, forwardChimeraWithK)
{
    real mfa = 1;
    real mfb = 1;
    real mfc = 1;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;

    vf::lbm::forwardChimeraWithK(mfa, mfb, mfc, vv, v2, K);

    EXPECT_THAT(mfa, REAL_EQ(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, REAL_EQ(-4.)); // -(mfa + mfb + mfc)
    EXPECT_THAT(mfc, REAL_EQ(6.));  // (mfa + mfc) + (mfa + mfb + mfc)
}


TEST(ChimeraTest, backwardChimeraWithK)
{
    // starting with the result values from the test above.
    real mfa = 3.;
    real mfb = -4.;
    real mfc = 6.;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;

    vf::lbm::backwardChimeraWithK(mfa, mfb, mfc, vv, v2, K);

    // resulting in the start values from the test above.
    EXPECT_THAT(mfa, REAL_EQ(1.));
    EXPECT_THAT(mfb, REAL_EQ(1.));
    EXPECT_THAT(mfc, REAL_EQ(1.));
}

//! \}
