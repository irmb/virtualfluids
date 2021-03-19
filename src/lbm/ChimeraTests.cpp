#include <gmock/gmock.h>

#include "Chimera.h"

/*
* InverseChimeraWithK
*/
TEST(ChimeraTest, forwardInverseChimeraWithK)
{
    real mfa = 1;
    real mfb = 1;
    real mfc = 1;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;
    const real Kinverse = 1 / K;

    VF::LBM::forwardInverseChimeraWithK(mfa, mfb, mfc, vv, v2, K, Kinverse);

    EXPECT_THAT(mfa, testing::DoubleEq(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, testing::DoubleEq(-4.)); // -(mfa + mfb + mfc + 1)
    EXPECT_THAT(mfc, testing::DoubleEq(6.));  // (mfa + mfc) + (mfa + mfb + mfc + 1)
}


TEST(ChimeraTest, backwardInverseChimeraWithK)
{
    // starting with the result values from the test above.
    real mfa = 3.;
    real mfb = -4.;
    real mfc = 6.;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;
    const real Kinverse = 1 / K;

    VF::LBM::backwardInverseChimeraWithK(mfa, mfb, mfc, vv, v2, K, Kinverse);

    // resulting in the start values from the test above.
    EXPECT_THAT(mfa, testing::DoubleEq(1.));
    EXPECT_THAT(mfb, testing::DoubleEq(1.));
    EXPECT_THAT(mfc, testing::DoubleEq(1.));
}

/*
* Chimera
*/
TEST(ChimeraTest, forwardChimera)
{
    real mfa = 1;
    real mfb = 1;
    real mfc = 1;

    const real vv = 1.;
    const real v2 = 1.;

    VF::LBM::forwardChimera(mfa, mfb, mfc, vv, v2);

    EXPECT_THAT(mfa, testing::DoubleEq(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, testing::DoubleEq(-3.)); // -(mfa + mfb + mfc)
    EXPECT_THAT(mfc, testing::DoubleEq(5.));  // (mfa + mfc) + (mfa + mfb + mfc)
}


TEST(ChimeraTest, backwardChimera)
{
    // starting with the result values from the test above.
    real mfa = 3.;
    real mfb = -3.;
    real mfc = 5.;

    const real vv = 1.;
    const real v2 = 1.;

    VF::LBM::backwardChimera(mfa, mfb, mfc, vv, v2);

    // resulting in the start values from the test above.
    EXPECT_THAT(mfa, testing::DoubleEq(1.));
    EXPECT_THAT(mfb, testing::DoubleEq(1.));
    EXPECT_THAT(mfc, testing::DoubleEq(1.));
}

/*
* ChimeraWithK
*/
TEST(ChimeraTest, forwardChimeraWithK)
{
    real mfa = 1;
    real mfb = 1;
    real mfc = 1;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;

    VF::LBM::forwardChimeraWithK(mfa, mfb, mfc, vv, v2, K);

    EXPECT_THAT(mfa, testing::DoubleEq(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, testing::DoubleEq(-4.)); // -(mfa + mfb + mfc)
    EXPECT_THAT(mfc, testing::DoubleEq(6.));  // (mfa + mfc) + (mfa + mfb + mfc)
}


TEST(ChimeraTest, backwardChimeraWithK)
{
    // starting with the result values from the test above.
    real mfa = 3.;
    real mfb = -4.;
    real mfc = 6.;

    const real vv = 1.;
    const real v2 = 1.;

    const real K = 1.;

    VF::LBM::backwardChimeraWithK(mfa, mfb, mfc, vv, v2, K);

    // resulting in the start values from the test above.
    EXPECT_THAT(mfa, testing::DoubleEq(1.));
    EXPECT_THAT(mfb, testing::DoubleEq(1.));
    EXPECT_THAT(mfc, testing::DoubleEq(1.));
}
