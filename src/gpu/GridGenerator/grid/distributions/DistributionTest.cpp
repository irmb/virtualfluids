# include <gmock/gmock.h>
#include <string>
# include "Distribution.h"

#include "grid/distributions/D3Q27.h"

TEST(DistributionTest, DistributionReturnsCorrectDirections)
{
    Distribution dist = DistributionHelper::getDistribution27();

    EXPECT_THAT(dist.directions[DIR_27_E][0], testing::Eq(DIR_27_E_X));
    EXPECT_THAT(dist.directions[DIR_27_E][1], testing::Eq(DIR_27_E_Y));
    EXPECT_THAT(dist.directions[DIR_27_E][2], testing::Eq(DIR_27_E_Z));
    EXPECT_THAT(dist.dirs[DIR_27_E * 3    ], testing::Eq(DIR_27_E_X));
    EXPECT_THAT(dist.dirs[DIR_27_E * 3 + 1], testing::Eq(DIR_27_E_Y));
    EXPECT_THAT(dist.dirs[DIR_27_E * 3 + 2], testing::Eq(DIR_27_E_Z));

    EXPECT_THAT(dist.directions[DIR_27_B][0], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.directions[DIR_27_B][1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.directions[DIR_27_B][2], testing::Eq(DIR_27_B_Z));
    EXPECT_THAT(dist.dirs[DIR_27_B * 3    ], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.dirs[DIR_27_B * 3 + 1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.dirs[DIR_27_B * 3 + 2], testing::Eq(DIR_27_B_Z));
    
    EXPECT_THAT(dist.directions[DIR_27_REST][0], testing::Eq(0));
    EXPECT_THAT(dist.directions[DIR_27_REST][1], testing::Eq(0));
    EXPECT_THAT(dist.directions[DIR_27_REST][2], testing::Eq(0));
    EXPECT_THAT(dist.dirs[DIR_27_REST * 3    ], testing::Eq(0));
    EXPECT_THAT(dist.dirs[DIR_27_REST * 3 + 1], testing::Eq(0));
    EXPECT_THAT(dist.dirs[DIR_27_REST * 3 + 2], testing::Eq(0));

    EXPECT_THAT(dist.directions[DIR_27_NE][0], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.directions[DIR_27_NE][1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.directions[DIR_27_NE][2], testing::Eq(DIR_27_NE_Z));
    EXPECT_THAT(dist.dirs[DIR_27_NE * 3    ], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.dirs[DIR_27_NE * 3 + 1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.dirs[DIR_27_NE * 3 + 2], testing::Eq(DIR_27_NE_Z));

    EXPECT_THAT(dist.directions[DIR_27_TS][0], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.directions[DIR_27_TS][1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.directions[DIR_27_TS][2], testing::Eq(DIR_27_TS_Z));
    EXPECT_THAT(dist.dirs[DIR_27_TS * 3    ], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.dirs[DIR_27_TS * 3 + 1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.dirs[DIR_27_TS * 3 + 2], testing::Eq(DIR_27_TS_Z));

    EXPECT_THAT(dist.directions[DIR_27_TNE][0], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.directions[DIR_27_TNE][1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.directions[DIR_27_TNE][2], testing::Eq(DIR_27_TNE_Z));
    EXPECT_THAT(dist.dirs[DIR_27_TNE * 3    ], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.dirs[DIR_27_TNE * 3 + 1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.dirs[DIR_27_TNE * 3 + 2], testing::Eq(DIR_27_TNE_Z));

    EXPECT_THAT(dist.directions[DIR_27_BSW][0], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.directions[DIR_27_BSW][1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.directions[DIR_27_BSW][2], testing::Eq(DIR_27_BSW_Z));
    EXPECT_THAT(dist.dirs[DIR_27_BSW * 3    ], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.dirs[DIR_27_BSW * 3 + 1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.dirs[DIR_27_BSW * 3 + 2], testing::Eq(DIR_27_BSW_Z));
}