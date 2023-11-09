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

    EXPECT_THAT(dist.directions[d00M][0], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.directions[d00M][1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.directions[d00M][2], testing::Eq(DIR_27_B_Z));
    EXPECT_THAT(dist.dirs[d00M * 3    ], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.dirs[d00M * 3 + 1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.dirs[d00M * 3 + 2], testing::Eq(DIR_27_B_Z));
    
    EXPECT_THAT(dist.directions[d000][0], testing::Eq(0));
    EXPECT_THAT(dist.directions[d000][1], testing::Eq(0));
    EXPECT_THAT(dist.directions[d000][2], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3    ], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3 + 1], testing::Eq(0));
    EXPECT_THAT(dist.dirs[d000 * 3 + 2], testing::Eq(0));

    EXPECT_THAT(dist.directions[dPP0][0], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.directions[dPP0][1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.directions[dPP0][2], testing::Eq(DIR_27_NE_Z));
    EXPECT_THAT(dist.dirs[dPP0 * 3    ], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.dirs[dPP0 * 3 + 1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.dirs[dPP0 * 3 + 2], testing::Eq(DIR_27_NE_Z));

    EXPECT_THAT(dist.directions[d0MP][0], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.directions[d0MP][1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.directions[d0MP][2], testing::Eq(DIR_27_TS_Z));
    EXPECT_THAT(dist.dirs[d0MP * 3    ], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.dirs[d0MP * 3 + 1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.dirs[d0MP * 3 + 2], testing::Eq(DIR_27_TS_Z));

    EXPECT_THAT(dist.directions[dPPP][0], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.directions[dPPP][1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.directions[dPPP][2], testing::Eq(DIR_27_TNE_Z));
    EXPECT_THAT(dist.dirs[dPPP * 3    ], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.dirs[dPPP * 3 + 1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.dirs[dPPP * 3 + 2], testing::Eq(DIR_27_TNE_Z));

    EXPECT_THAT(dist.directions[dMMM][0], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.directions[dMMM][1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.directions[dMMM][2], testing::Eq(DIR_27_BSW_Z));
    EXPECT_THAT(dist.dirs[dMMM * 3    ], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.dirs[dMMM * 3 + 1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.dirs[dMMM * 3 + 2], testing::Eq(DIR_27_BSW_Z));
}