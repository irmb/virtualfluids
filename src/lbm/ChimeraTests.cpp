#include <gmock/gmock.h>

#include "Chimera.h"
#include "D3Q27.h"



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
