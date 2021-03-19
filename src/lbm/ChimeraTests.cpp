#include <gmock/gmock.h>

#include "Chimera.h"

#ifdef VF_DOUBLE_ACCURACY
#define REAL_EQ(a) testing::DoubleEq(a)
#else
#define REAL_EQ(a) testing::FloatEq(a)
#endif

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

    EXPECT_THAT(mfa, REAL_EQ(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, REAL_EQ(-4.)); // -(mfa + mfb + mfc + 1)
    EXPECT_THAT(mfc, REAL_EQ(6.));  // (mfa + mfc) + (mfa + mfb + mfc + 1)
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
    EXPECT_THAT(mfa, REAL_EQ(1.));
    EXPECT_THAT(mfb, REAL_EQ(1.));
    EXPECT_THAT(mfc, REAL_EQ(1.));
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

    EXPECT_THAT(mfa, REAL_EQ(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, REAL_EQ(-3.)); // -(mfa + mfb + mfc)
    EXPECT_THAT(mfc, REAL_EQ(5.));  // (mfa + mfc) + (mfa + mfb + mfc)
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
    EXPECT_THAT(mfa, REAL_EQ(1.));
    EXPECT_THAT(mfb, REAL_EQ(1.));
    EXPECT_THAT(mfc, REAL_EQ(1.));
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

    EXPECT_THAT(mfa, REAL_EQ(3.));  // mfa + mfb + mfc
    EXPECT_THAT(mfb, REAL_EQ(-4.)); // -(mfa + mfb + mfc)
    EXPECT_THAT(mfc, REAL_EQ(6.));  // (mfa + mfc) + (mfa + mfb + mfc)
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
    EXPECT_THAT(mfa, REAL_EQ(1.));
    EXPECT_THAT(mfb, REAL_EQ(1.));
    EXPECT_THAT(mfc, REAL_EQ(1.));
}
