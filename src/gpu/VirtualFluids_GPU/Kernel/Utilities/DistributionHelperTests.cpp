#include <gmock/gmock.h>

#include "DistributionHelper.cuh"

#include "lbm/constants/D3Q27.h"
using namespace vf::lbm::dir;

auto RealEq = [](auto value) { 
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value); 
#else 
    return testing::FloatEq(value);
#endif
};


TEST(DistributionHelperTests, getPointerToDistribution_WhenEvenTimeStep_ShouldBeEqualToInput)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const uint size_Mat = 1;
    const bool isEvenTimeStep = true;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[E], RealEq(distributions_in[E]));
    EXPECT_THAT(*distribution_out.f[W], RealEq(distributions_in[W]));
    EXPECT_THAT(*distribution_out.f[N], RealEq(distributions_in[N]));
    EXPECT_THAT(*distribution_out.f[S], RealEq(distributions_in[S]));
    EXPECT_THAT(*distribution_out.f[T], RealEq(distributions_in[T]));
    EXPECT_THAT(*distribution_out.f[B], RealEq(distributions_in[B]));
    EXPECT_THAT(*distribution_out.f[NE], RealEq(distributions_in[NE]));
    EXPECT_THAT(*distribution_out.f[SW], RealEq(distributions_in[SW]));
    EXPECT_THAT(*distribution_out.f[SE], RealEq(distributions_in[SE]));
    EXPECT_THAT(*distribution_out.f[NW], RealEq(distributions_in[NW]));
    EXPECT_THAT(*distribution_out.f[TE], RealEq(distributions_in[TE]));
    EXPECT_THAT(*distribution_out.f[BW], RealEq(distributions_in[BW]));
    EXPECT_THAT(*distribution_out.f[BE], RealEq(distributions_in[BE]));
    EXPECT_THAT(*distribution_out.f[TW], RealEq(distributions_in[TW]));
    EXPECT_THAT(*distribution_out.f[TN], RealEq(distributions_in[TN]));
    EXPECT_THAT(*distribution_out.f[BS], RealEq(distributions_in[BS]));
    EXPECT_THAT(*distribution_out.f[BN], RealEq(distributions_in[BN]));
    EXPECT_THAT(*distribution_out.f[TS], RealEq(distributions_in[TS]));
    EXPECT_THAT(*distribution_out.f[REST], RealEq(distributions_in[REST]));
    EXPECT_THAT(*distribution_out.f[TNE], RealEq(distributions_in[TNE]));
    EXPECT_THAT(*distribution_out.f[TSW], RealEq(distributions_in[TSW]));
    EXPECT_THAT(*distribution_out.f[TSE], RealEq(distributions_in[TSE]));
    EXPECT_THAT(*distribution_out.f[TNW], RealEq(distributions_in[TNW]));
    EXPECT_THAT(*distribution_out.f[BNE], RealEq(distributions_in[BNE]));
    EXPECT_THAT(*distribution_out.f[BSW], RealEq(distributions_in[BSW]));
    EXPECT_THAT(*distribution_out.f[BSE], RealEq(distributions_in[BSE]));
    EXPECT_THAT(*distribution_out.f[BNW], RealEq(distributions_in[BNW]));
}

TEST(DistributionHelperTests, getPointerToDistribution_WhenOddTimeStep_ShouldBeSwapped)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const int size_Mat = 1;
    const bool isEvenTimeStep = false;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[W], RealEq(distributions_in[E]));
    EXPECT_THAT(*distribution_out.f[E], RealEq(distributions_in[W]));
    EXPECT_THAT(*distribution_out.f[S], RealEq(distributions_in[N]));
    EXPECT_THAT(*distribution_out.f[N], RealEq(distributions_in[S]));
    EXPECT_THAT(*distribution_out.f[B], RealEq(distributions_in[T]));
    EXPECT_THAT(*distribution_out.f[T], RealEq(distributions_in[B]));
    EXPECT_THAT(*distribution_out.f[SW], RealEq(distributions_in[NE]));
    EXPECT_THAT(*distribution_out.f[NE], RealEq(distributions_in[SW]));
    EXPECT_THAT(*distribution_out.f[NW], RealEq(distributions_in[SE]));
    EXPECT_THAT(*distribution_out.f[SE], RealEq(distributions_in[NW]));
    EXPECT_THAT(*distribution_out.f[BW], RealEq(distributions_in[TE]));
    EXPECT_THAT(*distribution_out.f[TE], RealEq(distributions_in[BW]));
    EXPECT_THAT(*distribution_out.f[TW], RealEq(distributions_in[BE]));
    EXPECT_THAT(*distribution_out.f[BE], RealEq(distributions_in[TW]));
    EXPECT_THAT(*distribution_out.f[BS], RealEq(distributions_in[TN]));
    EXPECT_THAT(*distribution_out.f[TN], RealEq(distributions_in[BS]));
    EXPECT_THAT(*distribution_out.f[TS], RealEq(distributions_in[BN]));
    EXPECT_THAT(*distribution_out.f[BN], RealEq(distributions_in[TS]));
    EXPECT_THAT(*distribution_out.f[REST], RealEq(distributions_in[REST]));
    EXPECT_THAT(*distribution_out.f[BSW], RealEq(distributions_in[TNE]));
    EXPECT_THAT(*distribution_out.f[BNE], RealEq(distributions_in[TSW]));
    EXPECT_THAT(*distribution_out.f[BNW], RealEq(distributions_in[TSE]));
    EXPECT_THAT(*distribution_out.f[BSE], RealEq(distributions_in[TNW]));
    EXPECT_THAT(*distribution_out.f[TSW], RealEq(distributions_in[BNE]));
    EXPECT_THAT(*distribution_out.f[TNE], RealEq(distributions_in[BSW]));
    EXPECT_THAT(*distribution_out.f[TNW], RealEq(distributions_in[BSE]));
    EXPECT_THAT(*distribution_out.f[TSE], RealEq(distributions_in[BNW]));
}
