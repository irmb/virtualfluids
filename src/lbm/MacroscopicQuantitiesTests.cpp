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
//! \author Soeren Peters
//=======================================================================================
#include <gmock/gmock.h>

#include <array>

#include "MacroscopicQuantities.h"
#include "constants/D3Q27.h"

#include <basics/tests/testUtilities.h>


/*
* given distributions, which are all one.
*/
// real f[27] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

class MacroscopicQuantitiesTest : public testing::Test {

    void SetUp() override {
       f = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    }

    public:
    std::array<real, 27> f;
};

using namespace vf::lbm;


TEST_F(MacroscopicQuantitiesTest, check_density)
{
    const double density = getDensity(f.data());

    const double expected_density = 27.;
    ASSERT_THAT(density, testing::DoubleEq(expected_density));
}

TEST_F(MacroscopicQuantitiesTest, whenFsAreEqual_velocityInEachDirectionShouldBeZero)
{
    const double velocityX1 = getIncompressibleVelocityX1(f.data());
    const double velocityX2 = getIncompressibleVelocityX2(f.data());
    const double velocityX3 = getIncompressibleVelocityX3(f.data());

    const double expected_velocity = 0.;
    EXPECT_THAT(velocityX1, testing::DoubleEq(expected_velocity));
    EXPECT_THAT(velocityX2, testing::DoubleEq(expected_velocity));
    EXPECT_THAT(velocityX3, testing::DoubleEq(expected_velocity));
}

TEST_F(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Eis2_velocityInX1ShouldBeOne)
{
    f[dir::dP00] = 2.;

    const double velocityX1 = getIncompressibleVelocityX1(f.data());
    const double velocityX2 = getIncompressibleVelocityX2(f.data());
    const double velocityX3 = getIncompressibleVelocityX3(f.data());

    const double expected_velocity_x1 = 1.;
    const double expected_velocity_x2 = 0.;
    const double expected_velocity_x3 = 0.;

    EXPECT_THAT(velocityX1, testing::DoubleEq(expected_velocity_x1));
    EXPECT_THAT(velocityX2, testing::DoubleEq(expected_velocity_x2));
    EXPECT_THAT(velocityX3, testing::DoubleEq(expected_velocity_x3));
}

TEST_F(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Nis2_velocityInX2ShouldBeOne)
{
    f[dir::d0P0] = 2.;

    const double velocity = getIncompressibleVelocityX2(f.data());

    const double expected_velocity = 1.;
    ASSERT_THAT(velocity, testing::DoubleEq(expected_velocity));
}


TEST_F(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Tis2_velocityInX3ShouldBeOne)
{
    f[dir::d00P] = 2.;

    const double velocity = getIncompressibleVelocityX3(f.data());

    const double expected_velocity = 1.;
    ASSERT_THAT(velocity, testing::DoubleEq(expected_velocity));
}

TEST_F(MacroscopicQuantitiesTest, givenAllFsAreOne_checkCompressibleValues)
{
    f[dir::dP00] = 2.;

    real rho{0}, vx1{0}, vx2{0}, vx3{0};
    getCompressibleMacroscopicValues(f.data(), rho, vx1, vx2, vx3);

    std::cout << getIncompressibleVelocityX1(f.data()) << std::endl;

    const real expected_rho = 28.;
    const real expected_velocity = (1. / (1 + expected_rho));

    EXPECT_THAT(rho, RealEq(expected_rho));
    EXPECT_THAT(vx1, RealEq(expected_velocity));
    EXPECT_THAT(vx2, RealEq(0.));
    EXPECT_THAT(vx3, RealEq(0.));
}
