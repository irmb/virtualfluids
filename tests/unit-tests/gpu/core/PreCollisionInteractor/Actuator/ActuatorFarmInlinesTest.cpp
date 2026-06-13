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
//! \addtogroup gpu_PreCollisionInteractor_tests PreCollisionInteractor
//! \ingroup gpu_core_tests
//! \{
#include <gmock/gmock.h>

#include <gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmInlines.h>

using namespace vf::gpu;

TEST(ActuatorFarmInlinesTest, calcPointIndexInBladeArrays)
{
    const uint numberOfPointsPerBlade = 4;
    const uint numberOfBlades = 3;

    // first node on first blade
    uint bladePoint = 0;
    uint blade = 0;
    uint turbine = 0;
    auto pointIndexInBladeArrays = calcPointIndexInBladeArrays(bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(0));
    pointIndexInBladeArrays = calcPointIndexInBladeArrays({turbine, blade, bladePoint}, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(0));

    // last node on first blade
    bladePoint = 3;
    blade = 0;
    turbine = 0;
    pointIndexInBladeArrays = calcPointIndexInBladeArrays(bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(3));
    pointIndexInBladeArrays = calcPointIndexInBladeArrays({turbine, blade, bladePoint}, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(3));

    // first node on third blade
    bladePoint = 0;
    blade = 2;
    turbine = 0;
    pointIndexInBladeArrays = calcPointIndexInBladeArrays(bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(8));
    pointIndexInBladeArrays = calcPointIndexInBladeArrays({turbine, blade, bladePoint}, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(8));

    // last node on third blade, also last node on first turbine
    bladePoint = 3;
    blade = 2;
    turbine = 0;
    pointIndexInBladeArrays = calcPointIndexInBladeArrays(bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(11));
    pointIndexInBladeArrays = calcPointIndexInBladeArrays({turbine, blade, bladePoint}, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(11));

    // first node on second turbine
    bladePoint = 0;
    blade = 0;
    turbine = 1;
    pointIndexInBladeArrays = calcPointIndexInBladeArrays(bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(12));
    pointIndexInBladeArrays = calcPointIndexInBladeArrays({turbine, blade, bladePoint}, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(12));

    // last node on second turbine
    bladePoint = 3;
    blade = 2;
    turbine = 1;
    pointIndexInBladeArrays = calcPointIndexInBladeArrays(bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(23));
    pointIndexInBladeArrays = calcPointIndexInBladeArrays({turbine, blade, bladePoint}, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(pointIndexInBladeArrays, testing::Eq(23));
}

TEST(ActuatorFarmInlinesTest, calcTurbineBladeAndBladePoint)
{
    const uint numberOfPointsPerBlade = 4;
    const uint numberOfBlades = 3;

    uint bladePoint;
    uint blade;
    uint turbine;

    TurbinePointIndex result;

    uint node = 0; // first node on first blade
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladePoint, testing::Eq(0));
    EXPECT_THAT(blade, testing::Eq(0));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladePoint(node, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladePoint, testing::Eq(0));
    EXPECT_THAT(result.blade, testing::Eq(0));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 3; // last node on first blade
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladePoint, testing::Eq(3));
    EXPECT_THAT(blade, testing::Eq(0));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladePoint(node, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladePoint, testing::Eq(3));
    EXPECT_THAT(result.blade, testing::Eq(0));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 8; // first node on third blade
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladePoint, testing::Eq(0));
    EXPECT_THAT(blade, testing::Eq(2));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladePoint(node, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladePoint, testing::Eq(0));
    EXPECT_THAT(result.blade, testing::Eq(2));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 11; // last node on third blade, also last node on first turbine
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladePoint, testing::Eq(3));
    EXPECT_THAT(blade, testing::Eq(2));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladePoint(node, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladePoint, testing::Eq(3));
    EXPECT_THAT(result.blade, testing::Eq(2));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 12; // first node on second turbine
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladePoint, testing::Eq(0));
    EXPECT_THAT(blade, testing::Eq(0));
    EXPECT_THAT(turbine, testing::Eq(1));
    result = calcTurbineBladeAndBladePoint(node, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladePoint, testing::Eq(0));
    EXPECT_THAT(result.blade, testing::Eq(0));
    EXPECT_THAT(result.turbine, testing::Eq(1));

    node = 23; // last node on second turbine
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladePoint, testing::Eq(3));
    EXPECT_THAT(blade, testing::Eq(2));
    EXPECT_THAT(turbine, testing::Eq(1));
    result = calcTurbineBladeAndBladePoint(node, numberOfPointsPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladePoint, testing::Eq(3));
    EXPECT_THAT(result.blade, testing::Eq(2));
    EXPECT_THAT(result.turbine, testing::Eq(1));
}

//! \}
