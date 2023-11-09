#include "interpolateStaticToRotatingInlines.h"
#include "tests/testUtilities.h"

class RotateSecondOrderMomentsFromRotatingToGlobalTest : public testing::Test
{
protected:
    const real mxxMyyBeforeRotation = -0.15;
    const real mxxMzzBeforeRotation = -0.21;
    const real m011BeforeRotation = 0.1;
    const real m101BeforeRotation = 0.22;
    const real m110BeforeRotation = 0.333;

    real mxxMyy = mxxMyyBeforeRotation;
    real mxxMzz = mxxMzzBeforeRotation;
    real m011 = m011BeforeRotation;
    real m101 = m101BeforeRotation;
    real m110 = m110BeforeRotation;
};

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, angleIsZero_momentsDoNotChange)
{
    std::array<real, 3> angles = { 0.0, 0.0, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    EXPECT_THAT(mxxMyy, RealEq(mxxMyyBeforeRotation));
    EXPECT_THAT(mxxMzz, RealEq(mxxMzzBeforeRotation));
    EXPECT_THAT(m011, RealEq(m011BeforeRotation));
    EXPECT_THAT(m101, RealEq(m101BeforeRotation));
    EXPECT_THAT(m110, RealEq(m110BeforeRotation));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, anglesAre2Pi_momentsDoNotChange)
{
    std::array<real, 3> angles = { c2Pi, c2Pi, c2Pi };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m011, RealNear(m011BeforeRotation, 1e-6));
    EXPECT_THAT(m101, RealNear(m101BeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, xAngleIsPi_someMomentsReverse)
{
    std::array<real, 3> angles = { cPi, 0.0, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // same
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m011, RealNear(m011BeforeRotation, 1e-6));

    // reversed
    EXPECT_THAT(m101, RealNear(-m101BeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(-m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, yAngleIsPi_someMomentsReverse)
{
    std::array<real, 3> angles = { 0.0, cPi, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // same
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m101, RealNear(m101BeforeRotation, 1e-6));

    // reversed
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(-m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, zAngleIsPi_someMomentsReverse)
{
    std::array<real, 3> angles = { 0.0, 0.0, cPi };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // same
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(m110BeforeRotation, 1e-6));

    // reversed
    EXPECT_THAT(m101, RealNear(-m101BeforeRotation, 1e-6));
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, xAngleIsPiHalf_momentsChange)
{
    std::array<real, 3> angles = { (real)0.5 * cPi, 0.0, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);


    // reversed
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(mxxMyy, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(m101, RealNear(m110BeforeRotation, 1e-6));
    
    // switched and reversed
    EXPECT_THAT(m110, RealNear(-m101BeforeRotation, 1e-6));

    real mxxPyyPzz = 0.0;
    real useNEQ = 1.0;
    real m200BeforeRotation = c1o3 * (        mxxMyyBeforeRotation +        mxxMzzBeforeRotation + mxxPyyPzz) * useNEQ;
    real m020BeforeRotation = c1o3 * (-c2o1 * mxxMyyBeforeRotation +        mxxMzzBeforeRotation + mxxPyyPzz) * useNEQ;
    real m002BeforeRotation = c1o3 * (        mxxMyyBeforeRotation - c2o1 * mxxMzzBeforeRotation + mxxPyyPzz) * useNEQ;
    real m200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    real m020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    real m002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    EXPECT_THAT(m200, RealNear(m200BeforeRotation, 1e-6));
    EXPECT_THAT(m020, RealNear(m002BeforeRotation, 1e-6));
    EXPECT_THAT(m002, RealNear(m020BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, xAngleIsMinusPiHalf_momentsChange)
{
    std::array<real, 3> angles = { (real)-0.5 * cPi, 0.0, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(mxxMyy, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(m101BeforeRotation, 1e-6));
    
    // switched and reversed
    EXPECT_THAT(m101, RealNear(-m110BeforeRotation, 1e-6));

    real mxxPyyPzz = 0.0;
    real useNEQ = 1.0;
    real m200BeforeRotation = c1o3 * (        mxxMyyBeforeRotation +        mxxMzzBeforeRotation + mxxPyyPzz) * useNEQ;
    real m020BeforeRotation = c1o3 * (-c2o1 * mxxMyyBeforeRotation +        mxxMzzBeforeRotation + mxxPyyPzz) * useNEQ;
    real m002BeforeRotation = c1o3 * (        mxxMyyBeforeRotation - c2o1 * mxxMzzBeforeRotation + mxxPyyPzz) * useNEQ;
    real m200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    real m020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    real m002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    EXPECT_THAT(m200, RealNear(m200BeforeRotation, 1e-6));
    EXPECT_THAT(m020, RealNear(m002BeforeRotation, 1e-6));
    EXPECT_THAT(m002, RealNear(m020BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, yAngleIsPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, (real)0.5 * cPi, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(m101, RealNear(-m101BeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(-mxxMzzBeforeRotation, 1e-6));

    // combined
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation - mxxMzzBeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(m110, RealNear(m011BeforeRotation, 1e-6));

    // switched and reversed
    EXPECT_THAT(m011, RealNear(-m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, yAngleIsMinusPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, (real)-0.5 * cPi, 0.0 };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(m101, RealNear(-m101BeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(-mxxMzzBeforeRotation, 1e-6));

    // combined
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation - mxxMzzBeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(m011, RealNear(m110BeforeRotation, 1e-6));

    // switched and reversed
    EXPECT_THAT(m110, RealNear(-m011BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, zAngleIsPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, 0.0, (real)0.5 * cPi};

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(m110, RealNear(-m110BeforeRotation, 1e-6));
    EXPECT_THAT(mxxMyy, RealNear(-mxxMyyBeforeRotation, 1e-6));

    // combined
    EXPECT_THAT(mxxMzz, RealNear(-mxxMyyBeforeRotation + mxxMzzBeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(m011, RealNear(m101BeforeRotation, 1e-6));

    // switched and reversed
    EXPECT_THAT(m101, RealNear(-m011BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromRotatingToGlobalTest, zAngleIsMinusPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, 0.0, (real)-0.5 * cPi };

    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(m110, RealNear(-m110BeforeRotation, 1e-6));
    EXPECT_THAT(mxxMyy, RealNear(-mxxMyyBeforeRotation, 1e-6));

    // combined
    EXPECT_THAT(mxxMzz, RealNear(-mxxMyyBeforeRotation + mxxMzzBeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(m101, RealNear(m011BeforeRotation, 1e-6));

    // switched and reversed
    EXPECT_THAT(m011, RealNear(-m101BeforeRotation, 1e-6));
}
