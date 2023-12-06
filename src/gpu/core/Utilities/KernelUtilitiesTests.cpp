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
//! \author Martin Schoenherr
//=======================================================================================
#include <gmock/gmock.h>

#include <basics/tests/testUtilities.h>

#include "KernelUtilities.h"

#include <lbm/constants/D3Q27.h>

using namespace vf::lbm::dir;

TEST(DistributionHelperTests, getPointerToDistribution_WhenEvenTimeStep_ShouldBeEqualToInput)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const uint size_Mat = 1;
    const bool isEvenTimeStep = true;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[dP00], RealEq(distributions_in[dP00]));
    EXPECT_THAT(*distribution_out.f[dM00], RealEq(distributions_in[dM00]));
    EXPECT_THAT(*distribution_out.f[d0P0], RealEq(distributions_in[d0P0]));
    EXPECT_THAT(*distribution_out.f[d0M0], RealEq(distributions_in[d0M0]));
    EXPECT_THAT(*distribution_out.f[d00P], RealEq(distributions_in[d00P]));
    EXPECT_THAT(*distribution_out.f[d00M], RealEq(distributions_in[d00M]));
    EXPECT_THAT(*distribution_out.f[dPP0], RealEq(distributions_in[dPP0]));
    EXPECT_THAT(*distribution_out.f[dMM0], RealEq(distributions_in[dMM0]));
    EXPECT_THAT(*distribution_out.f[dPM0], RealEq(distributions_in[dPM0]));
    EXPECT_THAT(*distribution_out.f[dMP0], RealEq(distributions_in[dMP0]));
    EXPECT_THAT(*distribution_out.f[dP0P], RealEq(distributions_in[dP0P]));
    EXPECT_THAT(*distribution_out.f[dM0M], RealEq(distributions_in[dM0M]));
    EXPECT_THAT(*distribution_out.f[dP0M], RealEq(distributions_in[dP0M]));
    EXPECT_THAT(*distribution_out.f[dM0P], RealEq(distributions_in[dM0P]));
    EXPECT_THAT(*distribution_out.f[d0PP], RealEq(distributions_in[d0PP]));
    EXPECT_THAT(*distribution_out.f[d0MM], RealEq(distributions_in[d0MM]));
    EXPECT_THAT(*distribution_out.f[d0PM], RealEq(distributions_in[d0PM]));
    EXPECT_THAT(*distribution_out.f[d0MP], RealEq(distributions_in[d0MP]));
    EXPECT_THAT(*distribution_out.f[d000], RealEq(distributions_in[d000]));
    EXPECT_THAT(*distribution_out.f[dPPP], RealEq(distributions_in[dPPP]));
    EXPECT_THAT(*distribution_out.f[dMMP], RealEq(distributions_in[dMMP]));
    EXPECT_THAT(*distribution_out.f[dPMP], RealEq(distributions_in[dPMP]));
    EXPECT_THAT(*distribution_out.f[dMPP], RealEq(distributions_in[dMPP]));
    EXPECT_THAT(*distribution_out.f[dPPM], RealEq(distributions_in[dPPM]));
    EXPECT_THAT(*distribution_out.f[dMMM], RealEq(distributions_in[dMMM]));
    EXPECT_THAT(*distribution_out.f[dPMM], RealEq(distributions_in[dPMM]));
    EXPECT_THAT(*distribution_out.f[dMPM], RealEq(distributions_in[dMPM]));
}

TEST(DistributionHelperTests, getPointerToDistribution_WhenOddTimeStep_ShouldBeSwapped)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const int size_Mat = 1;
    const bool isEvenTimeStep = false;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[dM00], RealEq(distributions_in[dP00]));
    EXPECT_THAT(*distribution_out.f[dP00], RealEq(distributions_in[dM00]));
    EXPECT_THAT(*distribution_out.f[d0M0], RealEq(distributions_in[d0P0]));
    EXPECT_THAT(*distribution_out.f[d0P0], RealEq(distributions_in[d0M0]));
    EXPECT_THAT(*distribution_out.f[d00M], RealEq(distributions_in[d00P]));
    EXPECT_THAT(*distribution_out.f[d00P], RealEq(distributions_in[d00M]));
    EXPECT_THAT(*distribution_out.f[dMM0], RealEq(distributions_in[dPP0]));
    EXPECT_THAT(*distribution_out.f[dPP0], RealEq(distributions_in[dMM0]));
    EXPECT_THAT(*distribution_out.f[dMP0], RealEq(distributions_in[dPM0]));
    EXPECT_THAT(*distribution_out.f[dPM0], RealEq(distributions_in[dMP0]));
    EXPECT_THAT(*distribution_out.f[dM0M], RealEq(distributions_in[dP0P]));
    EXPECT_THAT(*distribution_out.f[dP0P], RealEq(distributions_in[dM0M]));
    EXPECT_THAT(*distribution_out.f[dM0P], RealEq(distributions_in[dP0M]));
    EXPECT_THAT(*distribution_out.f[dP0M], RealEq(distributions_in[dM0P]));
    EXPECT_THAT(*distribution_out.f[d0MM], RealEq(distributions_in[d0PP]));
    EXPECT_THAT(*distribution_out.f[d0PP], RealEq(distributions_in[d0MM]));
    EXPECT_THAT(*distribution_out.f[d0MP], RealEq(distributions_in[d0PM]));
    EXPECT_THAT(*distribution_out.f[d0PM], RealEq(distributions_in[d0MP]));
    EXPECT_THAT(*distribution_out.f[d000], RealEq(distributions_in[d000]));
    EXPECT_THAT(*distribution_out.f[dMMM], RealEq(distributions_in[dPPP]));
    EXPECT_THAT(*distribution_out.f[dPPM], RealEq(distributions_in[dMMP]));
    EXPECT_THAT(*distribution_out.f[dMPM], RealEq(distributions_in[dPMP]));
    EXPECT_THAT(*distribution_out.f[dPMM], RealEq(distributions_in[dMPP]));
    EXPECT_THAT(*distribution_out.f[dMMP], RealEq(distributions_in[dPPM]));
    EXPECT_THAT(*distribution_out.f[dPPP], RealEq(distributions_in[dMMM]));
    EXPECT_THAT(*distribution_out.f[dMPP], RealEq(distributions_in[dPMM]));
    EXPECT_THAT(*distribution_out.f[dPMP], RealEq(distributions_in[dMPM]));
}
