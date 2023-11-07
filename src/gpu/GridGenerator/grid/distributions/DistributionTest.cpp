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

    EXPECT_THAT(dist.directions[DIR_00M][0], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.directions[DIR_00M][1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.directions[DIR_00M][2], testing::Eq(DIR_27_B_Z));
    EXPECT_THAT(dist.dirs[DIR_00M * 3    ], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.dirs[DIR_00M * 3 + 1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.dirs[DIR_00M * 3 + 2], testing::Eq(DIR_27_B_Z));
    
    EXPECT_THAT(dist.directions[d000][0], testing::Eq(0));
    EXPECT_THAT(dist.directions[d000][1], testing::Eq(0));
    EXPECT_THAT(dist.directions[d000][2], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3    ], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3 + 1], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3 + 2], testing::Eq(0));

    EXPECT_THAT(dist.directions[DIR_PP0][0], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.directions[DIR_PP0][1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.directions[DIR_PP0][2], testing::Eq(DIR_27_NE_Z));
    EXPECT_THAT(dist.dirs[DIR_PP0 * 3    ], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.dirs[DIR_PP0 * 3 + 1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.dirs[DIR_PP0 * 3 + 2], testing::Eq(DIR_27_NE_Z));

    EXPECT_THAT(dist.directions[DIR_0MP][0], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.directions[DIR_0MP][1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.directions[DIR_0MP][2], testing::Eq(DIR_27_TS_Z));
    EXPECT_THAT(dist.dirs[DIR_0MP * 3    ], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.dirs[DIR_0MP * 3 + 1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.dirs[DIR_0MP * 3 + 2], testing::Eq(DIR_27_TS_Z));

    EXPECT_THAT(dist.directions[DIR_PPP][0], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.directions[DIR_PPP][1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.directions[DIR_PPP][2], testing::Eq(DIR_27_TNE_Z));
    EXPECT_THAT(dist.dirs[DIR_PPP * 3    ], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.dirs[DIR_PPP * 3 + 1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.dirs[DIR_PPP * 3 + 2], testing::Eq(DIR_27_TNE_Z));

    EXPECT_THAT(dist.directions[DIR_MMM][0], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.directions[DIR_MMM][1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.directions[DIR_MMM][2], testing::Eq(DIR_27_BSW_Z));
    EXPECT_THAT(dist.dirs[DIR_MMM * 3    ], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.dirs[DIR_MMM * 3 + 1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.dirs[DIR_MMM * 3 + 2], testing::Eq(DIR_27_BSW_Z));
}