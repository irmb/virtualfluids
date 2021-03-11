#include <gmock/gmock.h>

#include "CalcMac.h"

TEST(CalcMacTest, calcDensity)
{
    real f[27] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    double density = LBM::getDensity(f);

    ASSERT_THAT(density, testing::DoubleEq(27));
}