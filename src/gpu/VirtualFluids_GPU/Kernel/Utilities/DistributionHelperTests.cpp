#include <gmock/gmock.h>

#include "DistributionHelper.cuh"

#include "LBM/D3Q27.h"


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

    EXPECT_THAT(*distribution_out.f[dirE], RealEq(distributions_in[dirE]));
    EXPECT_THAT(*distribution_out.f[dirW], RealEq(distributions_in[dirW]));
    EXPECT_THAT(*distribution_out.f[dirN], RealEq(distributions_in[dirN]));
    EXPECT_THAT(*distribution_out.f[dirS], RealEq(distributions_in[dirS]));
    EXPECT_THAT(*distribution_out.f[dirT], RealEq(distributions_in[dirT]));
    EXPECT_THAT(*distribution_out.f[dirB], RealEq(distributions_in[dirB]));
    EXPECT_THAT(*distribution_out.f[dirNE], RealEq(distributions_in[dirNE]));
    EXPECT_THAT(*distribution_out.f[dirSW], RealEq(distributions_in[dirSW]));
    EXPECT_THAT(*distribution_out.f[dirSE], RealEq(distributions_in[dirSE]));
    EXPECT_THAT(*distribution_out.f[dirNW], RealEq(distributions_in[dirNW]));
    EXPECT_THAT(*distribution_out.f[dirTE], RealEq(distributions_in[dirTE]));
    EXPECT_THAT(*distribution_out.f[dirBW], RealEq(distributions_in[dirBW]));
    EXPECT_THAT(*distribution_out.f[dirBE], RealEq(distributions_in[dirBE]));
    EXPECT_THAT(*distribution_out.f[dirTW], RealEq(distributions_in[dirTW]));
    EXPECT_THAT(*distribution_out.f[dirTN], RealEq(distributions_in[dirTN]));
    EXPECT_THAT(*distribution_out.f[dirBS], RealEq(distributions_in[dirBS]));
    EXPECT_THAT(*distribution_out.f[dirBN], RealEq(distributions_in[dirBN]));
    EXPECT_THAT(*distribution_out.f[dirTS], RealEq(distributions_in[dirTS]));
    EXPECT_THAT(*distribution_out.f[dirREST], RealEq(distributions_in[dirREST]));
    EXPECT_THAT(*distribution_out.f[dirTNE], RealEq(distributions_in[dirTNE]));
    EXPECT_THAT(*distribution_out.f[dirTSW], RealEq(distributions_in[dirTSW]));
    EXPECT_THAT(*distribution_out.f[dirTSE], RealEq(distributions_in[dirTSE]));
    EXPECT_THAT(*distribution_out.f[dirTNW], RealEq(distributions_in[dirTNW]));
    EXPECT_THAT(*distribution_out.f[dirBNE], RealEq(distributions_in[dirBNE]));
    EXPECT_THAT(*distribution_out.f[dirBSW], RealEq(distributions_in[dirBSW]));
    EXPECT_THAT(*distribution_out.f[dirBSE], RealEq(distributions_in[dirBSE]));
    EXPECT_THAT(*distribution_out.f[dirBNW], RealEq(distributions_in[dirBNW]));
}

TEST(DistributionHelperTests, getPointerToDistribution_WhenOddTimeStep_ShouldBeSwapped)
{
    real distributions_in[27];
    for (int i = 0; i < 27; i++)
        distributions_in[i] = i;
    const int size_Mat = 1;
    const bool isEvenTimeStep = false;

    Distributions27 distribution_out = vf::gpu::getDistributionReferences27(distributions_in, size_Mat, isEvenTimeStep);

    EXPECT_THAT(*distribution_out.f[dirW], RealEq(distributions_in[dirE]));
    EXPECT_THAT(*distribution_out.f[dirE], RealEq(distributions_in[dirW]));
    EXPECT_THAT(*distribution_out.f[dirS], RealEq(distributions_in[dirN]));
    EXPECT_THAT(*distribution_out.f[dirN], RealEq(distributions_in[dirS]));
    EXPECT_THAT(*distribution_out.f[dirB], RealEq(distributions_in[dirT]));
    EXPECT_THAT(*distribution_out.f[dirT], RealEq(distributions_in[dirB]));
    EXPECT_THAT(*distribution_out.f[dirSW], RealEq(distributions_in[dirNE]));
    EXPECT_THAT(*distribution_out.f[dirNE], RealEq(distributions_in[dirSW]));
    EXPECT_THAT(*distribution_out.f[dirNW], RealEq(distributions_in[dirSE]));
    EXPECT_THAT(*distribution_out.f[dirSE], RealEq(distributions_in[dirNW]));
    EXPECT_THAT(*distribution_out.f[dirBW], RealEq(distributions_in[dirTE]));
    EXPECT_THAT(*distribution_out.f[dirTE], RealEq(distributions_in[dirBW]));
    EXPECT_THAT(*distribution_out.f[dirTW], RealEq(distributions_in[dirBE]));
    EXPECT_THAT(*distribution_out.f[dirBE], RealEq(distributions_in[dirTW]));
    EXPECT_THAT(*distribution_out.f[dirBS], RealEq(distributions_in[dirTN]));
    EXPECT_THAT(*distribution_out.f[dirTN], RealEq(distributions_in[dirBS]));
    EXPECT_THAT(*distribution_out.f[dirTS], RealEq(distributions_in[dirBN]));
    EXPECT_THAT(*distribution_out.f[dirBN], RealEq(distributions_in[dirTS]));
    EXPECT_THAT(*distribution_out.f[dirREST], RealEq(distributions_in[dirREST]));
    EXPECT_THAT(*distribution_out.f[dirBSW], RealEq(distributions_in[dirTNE]));
    EXPECT_THAT(*distribution_out.f[dirBNE], RealEq(distributions_in[dirTSW]));
    EXPECT_THAT(*distribution_out.f[dirBNW], RealEq(distributions_in[dirTSE]));
    EXPECT_THAT(*distribution_out.f[dirBSE], RealEq(distributions_in[dirTNW]));
    EXPECT_THAT(*distribution_out.f[dirTSW], RealEq(distributions_in[dirBNE]));
    EXPECT_THAT(*distribution_out.f[dirTNE], RealEq(distributions_in[dirBSW]));
    EXPECT_THAT(*distribution_out.f[dirTNW], RealEq(distributions_in[dirBSE]));
    EXPECT_THAT(*distribution_out.f[dirTSE], RealEq(distributions_in[dirBNW]));
}
