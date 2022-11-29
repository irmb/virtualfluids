#include <gmock/gmock.h>
#include "../utilities/testUtilities.h"

#include "DistributionHelper.cuh"

#include "lbm/constants/D3Q27.h"
using namespace vf::lbm::dir;

TEST(DistributionHelperTests, getPointerToDistribution_WhenEvenTimeStep_ShouldBeEqualToInput)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const uint size_Mat = 1;
    const bool isEvenTimeStep = true;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[DIR_P00], RealEq(distributions_in[DIR_P00]));
    EXPECT_THAT(*distribution_out.f[DIR_M00], RealEq(distributions_in[DIR_M00]));
    EXPECT_THAT(*distribution_out.f[DIR_0P0], RealEq(distributions_in[DIR_0P0]));
    EXPECT_THAT(*distribution_out.f[DIR_0M0], RealEq(distributions_in[DIR_0M0]));
    EXPECT_THAT(*distribution_out.f[DIR_00P], RealEq(distributions_in[DIR_00P]));
    EXPECT_THAT(*distribution_out.f[DIR_00M], RealEq(distributions_in[DIR_00M]));
    EXPECT_THAT(*distribution_out.f[DIR_PP0], RealEq(distributions_in[DIR_PP0]));
    EXPECT_THAT(*distribution_out.f[DIR_MM0], RealEq(distributions_in[DIR_MM0]));
    EXPECT_THAT(*distribution_out.f[DIR_PM0], RealEq(distributions_in[DIR_PM0]));
    EXPECT_THAT(*distribution_out.f[DIR_MP0], RealEq(distributions_in[DIR_MP0]));
    EXPECT_THAT(*distribution_out.f[DIR_P0P], RealEq(distributions_in[DIR_P0P]));
    EXPECT_THAT(*distribution_out.f[DIR_M0M], RealEq(distributions_in[DIR_M0M]));
    EXPECT_THAT(*distribution_out.f[DIR_P0M], RealEq(distributions_in[DIR_P0M]));
    EXPECT_THAT(*distribution_out.f[DIR_M0P], RealEq(distributions_in[DIR_M0P]));
    EXPECT_THAT(*distribution_out.f[DIR_0PP], RealEq(distributions_in[DIR_0PP]));
    EXPECT_THAT(*distribution_out.f[DIR_0MM], RealEq(distributions_in[DIR_0MM]));
    EXPECT_THAT(*distribution_out.f[DIR_0PM], RealEq(distributions_in[DIR_0PM]));
    EXPECT_THAT(*distribution_out.f[DIR_0MP], RealEq(distributions_in[DIR_0MP]));
    EXPECT_THAT(*distribution_out.f[DIR_000], RealEq(distributions_in[DIR_000]));
    EXPECT_THAT(*distribution_out.f[DIR_PPP], RealEq(distributions_in[DIR_PPP]));
    EXPECT_THAT(*distribution_out.f[DIR_MMP], RealEq(distributions_in[DIR_MMP]));
    EXPECT_THAT(*distribution_out.f[DIR_PMP], RealEq(distributions_in[DIR_PMP]));
    EXPECT_THAT(*distribution_out.f[DIR_MPP], RealEq(distributions_in[DIR_MPP]));
    EXPECT_THAT(*distribution_out.f[DIR_PPM], RealEq(distributions_in[DIR_PPM]));
    EXPECT_THAT(*distribution_out.f[DIR_MMM], RealEq(distributions_in[DIR_MMM]));
    EXPECT_THAT(*distribution_out.f[DIR_PMM], RealEq(distributions_in[DIR_PMM]));
    EXPECT_THAT(*distribution_out.f[DIR_MPM], RealEq(distributions_in[DIR_MPM]));
}

TEST(DistributionHelperTests, getPointerToDistribution_WhenOddTimeStep_ShouldBeSwapped)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const int size_Mat = 1;
    const bool isEvenTimeStep = false;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[DIR_M00], RealEq(distributions_in[DIR_P00]));
    EXPECT_THAT(*distribution_out.f[DIR_P00], RealEq(distributions_in[DIR_M00]));
    EXPECT_THAT(*distribution_out.f[DIR_0M0], RealEq(distributions_in[DIR_0P0]));
    EXPECT_THAT(*distribution_out.f[DIR_0P0], RealEq(distributions_in[DIR_0M0]));
    EXPECT_THAT(*distribution_out.f[DIR_00M], RealEq(distributions_in[DIR_00P]));
    EXPECT_THAT(*distribution_out.f[DIR_00P], RealEq(distributions_in[DIR_00M]));
    EXPECT_THAT(*distribution_out.f[DIR_MM0], RealEq(distributions_in[DIR_PP0]));
    EXPECT_THAT(*distribution_out.f[DIR_PP0], RealEq(distributions_in[DIR_MM0]));
    EXPECT_THAT(*distribution_out.f[DIR_MP0], RealEq(distributions_in[DIR_PM0]));
    EXPECT_THAT(*distribution_out.f[DIR_PM0], RealEq(distributions_in[DIR_MP0]));
    EXPECT_THAT(*distribution_out.f[DIR_M0M], RealEq(distributions_in[DIR_P0P]));
    EXPECT_THAT(*distribution_out.f[DIR_P0P], RealEq(distributions_in[DIR_M0M]));
    EXPECT_THAT(*distribution_out.f[DIR_M0P], RealEq(distributions_in[DIR_P0M]));
    EXPECT_THAT(*distribution_out.f[DIR_P0M], RealEq(distributions_in[DIR_M0P]));
    EXPECT_THAT(*distribution_out.f[DIR_0MM], RealEq(distributions_in[DIR_0PP]));
    EXPECT_THAT(*distribution_out.f[DIR_0PP], RealEq(distributions_in[DIR_0MM]));
    EXPECT_THAT(*distribution_out.f[DIR_0MP], RealEq(distributions_in[DIR_0PM]));
    EXPECT_THAT(*distribution_out.f[DIR_0PM], RealEq(distributions_in[DIR_0MP]));
    EXPECT_THAT(*distribution_out.f[DIR_000], RealEq(distributions_in[DIR_000]));
    EXPECT_THAT(*distribution_out.f[DIR_MMM], RealEq(distributions_in[DIR_PPP]));
    EXPECT_THAT(*distribution_out.f[DIR_PPM], RealEq(distributions_in[DIR_MMP]));
    EXPECT_THAT(*distribution_out.f[DIR_MPM], RealEq(distributions_in[DIR_PMP]));
    EXPECT_THAT(*distribution_out.f[DIR_PMM], RealEq(distributions_in[DIR_MPP]));
    EXPECT_THAT(*distribution_out.f[DIR_MMP], RealEq(distributions_in[DIR_PPM]));
    EXPECT_THAT(*distribution_out.f[DIR_PPP], RealEq(distributions_in[DIR_MMM]));
    EXPECT_THAT(*distribution_out.f[DIR_MPP], RealEq(distributions_in[DIR_PMM]));
    EXPECT_THAT(*distribution_out.f[DIR_PMP], RealEq(distributions_in[DIR_MPM]));
}
