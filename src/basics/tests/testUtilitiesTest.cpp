#include "testUtilities.h"

TEST(isEqualWithAccuracy, test)
{
    const real accuracy = 1.0;
    const real expected = 0.0;

    EXPECT_TRUE(isEqualWithAccuracy( 0.0, expected, accuracy));
    EXPECT_TRUE(isEqualWithAccuracy( 0.999999, expected, accuracy));
    EXPECT_TRUE(isEqualWithAccuracy( -0.999999, expected, accuracy));
    EXPECT_FALSE(isEqualWithAccuracy( 1.000001, expected, accuracy));
    EXPECT_FALSE(isEqualWithAccuracy( -1.000001, expected, accuracy));
}
