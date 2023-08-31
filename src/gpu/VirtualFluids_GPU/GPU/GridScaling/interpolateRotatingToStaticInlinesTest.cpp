#include "interpolateRotatingToStaticInlines.h"
#include "tests/testUtilities.h"

class RotateSecondOrderMomentsFromGlobalToRotatingTest : public testing::Test
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

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, angleIsZero_momentsDoNotChange)
{
    std::array<real, 3> angles = { 0.0, 0.0, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    EXPECT_THAT(mxxMyy, RealEq(mxxMyyBeforeRotation));
    EXPECT_THAT(mxxMzz, RealEq(mxxMzzBeforeRotation));
    EXPECT_THAT(m011, RealEq(m011BeforeRotation));
    EXPECT_THAT(m101, RealEq(m101BeforeRotation));
    EXPECT_THAT(m110, RealEq(m110BeforeRotation));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, angleIs2Pi_momentsDoNotChange)
{
    std::array<real, 3> angles = { c2Pi, c2Pi, c2Pi };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m011, RealNear(m011BeforeRotation, 1e-6));
    EXPECT_THAT(m101, RealNear(m101BeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, xAngleIsPi_someMomentsReverse)
{
    std::array<real, 3> angles = { cPi, 0.0, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // same
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m011, RealNear(m011BeforeRotation, 1e-6));

    // reversed
    EXPECT_THAT(m101, RealNear(-m101BeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(-m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, yAngleIsPi_someMomentsReverse)
{
    std::array<real, 3> angles = { 0.0, cPi, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // same
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m101, RealNear(m101BeforeRotation, 1e-6));

    // reversed
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(-m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, zAngleIsPi_someMomentsReverse)
{
    std::array<real, 3> angles = { 0.0, 0.0, cPi };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // same
    EXPECT_THAT(mxxMyy, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(m110BeforeRotation, 1e-6));

    // reversed
    EXPECT_THAT(m101, RealNear(-m101BeforeRotation, 1e-6));
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, xAngleIsPiHalf_momentsChange)
{
    std::array<real, 3> angles = { (real)0.5 * cPi, 0.0, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(m011, RealNear(-m011BeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(mxxMyy, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(m110, RealNear(m101BeforeRotation, 1e-6));
    
    // switched and reversed
    EXPECT_THAT(m101, RealNear(-m110BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, xAngleIsMinusPiHalf_momentsChange)
{
    std::array<real, 3> angles = { (real)-0.5 * cPi, 0.0, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

    // reversed
    EXPECT_THAT(-m011, RealNear(m011BeforeRotation, 1e-6));

    // switched
    EXPECT_THAT(mxxMyy, RealNear(mxxMzzBeforeRotation, 1e-6));
    EXPECT_THAT(mxxMzz, RealNear(mxxMyyBeforeRotation, 1e-6));
    EXPECT_THAT(m101, RealNear(m110BeforeRotation, 1e-6));
    
    // switched and reversed
    EXPECT_THAT(m110, RealNear(-m101BeforeRotation, 1e-6));
}

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, yAngleIsPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, (real)0.5 * cPi, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

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

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, yAngleIsMinusPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, (real)-0.5 * cPi, 0.0 };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

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

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, zAngleIsPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, 0.0, (real)0.5 * cPi};

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

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

TEST_F(RotateSecondOrderMomentsFromGlobalToRotatingTest, zAngleIsMinusPiHalf_momentsChange)
{
    std::array<real, 3> angles = { 0.0, 0.0, (real)-0.5 * cPi };

    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angles[0], angles[1], angles[2]);

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
