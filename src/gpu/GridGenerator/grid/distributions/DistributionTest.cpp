# include <gmock/gmock.h>
#include <string>
# include "Distribution.h"

#include "grid/distributions/D3Q27.h"

TEST(DistributionTest, DistributionReturnsCorrectDirections)
{
    Distribution dist = DistributionHelper::getDistribution27();

    EXPECT_THAT(dist.directions[E][0], testing::Eq(DIR_27_E_X));
    EXPECT_THAT(dist.directions[E][1], testing::Eq(DIR_27_E_Y));
    EXPECT_THAT(dist.directions[E][2], testing::Eq(DIR_27_E_Z));
    EXPECT_THAT(dist.dirs[E * 3    ], testing::Eq(DIR_27_E_X));
    EXPECT_THAT(dist.dirs[E * 3 + 1], testing::Eq(DIR_27_E_Y));
    EXPECT_THAT(dist.dirs[E * 3 + 2], testing::Eq(DIR_27_E_Z));

    EXPECT_THAT(dist.directions[B][0], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.directions[B][1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.directions[B][2], testing::Eq(DIR_27_B_Z));
    EXPECT_THAT(dist.dirs[B * 3    ], testing::Eq(DIR_27_B_X));
    EXPECT_THAT(dist.dirs[B * 3 + 1], testing::Eq(DIR_27_B_Y));
    EXPECT_THAT(dist.dirs[B * 3 + 2], testing::Eq(DIR_27_B_Z));
    
    EXPECT_THAT(dist.directions[REST][0], testing::Eq(0));
    EXPECT_THAT(dist.directions[REST][1], testing::Eq(0));
    EXPECT_THAT(dist.directions[REST][2], testing::Eq(0));
    EXPECT_THAT(dist.dirs[REST * 3    ], testing::Eq(0));
    EXPECT_THAT(dist.dirs[REST * 3 + 1], testing::Eq(0));
    EXPECT_THAT(dist.dirs[REST * 3 + 2], testing::Eq(0));

    EXPECT_THAT(dist.directions[NE][0], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.directions[NE][1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.directions[NE][2], testing::Eq(DIR_27_NE_Z));
    EXPECT_THAT(dist.dirs[NE * 3    ], testing::Eq(DIR_27_NE_X));
    EXPECT_THAT(dist.dirs[NE * 3 + 1], testing::Eq(DIR_27_NE_Y));
    EXPECT_THAT(dist.dirs[NE * 3 + 2], testing::Eq(DIR_27_NE_Z));

    EXPECT_THAT(dist.directions[TS][0], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.directions[TS][1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.directions[TS][2], testing::Eq(DIR_27_TS_Z));
    EXPECT_THAT(dist.dirs[TS * 3    ], testing::Eq(DIR_27_TS_X));
    EXPECT_THAT(dist.dirs[TS * 3 + 1], testing::Eq(DIR_27_TS_Y));
    EXPECT_THAT(dist.dirs[TS * 3 + 2], testing::Eq(DIR_27_TS_Z));

    EXPECT_THAT(dist.directions[TNE][0], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.directions[TNE][1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.directions[TNE][2], testing::Eq(DIR_27_TNE_Z));
    EXPECT_THAT(dist.dirs[TNE * 3    ], testing::Eq(DIR_27_TNE_X));
    EXPECT_THAT(dist.dirs[TNE * 3 + 1], testing::Eq(DIR_27_TNE_Y));
    EXPECT_THAT(dist.dirs[TNE * 3 + 2], testing::Eq(DIR_27_TNE_Z));

    EXPECT_THAT(dist.directions[BSW][0], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.directions[BSW][1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.directions[BSW][2], testing::Eq(DIR_27_BSW_Z));
    EXPECT_THAT(dist.dirs[BSW * 3    ], testing::Eq(DIR_27_BSW_X));
    EXPECT_THAT(dist.dirs[BSW * 3 + 1], testing::Eq(DIR_27_BSW_Y));
    EXPECT_THAT(dist.dirs[BSW * 3 + 2], testing::Eq(DIR_27_BSW_Z));
}