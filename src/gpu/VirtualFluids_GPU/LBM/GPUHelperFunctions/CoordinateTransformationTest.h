
#include "CoordinateTransformation.h"
#include "tests/testUtilities.h"

bool isEqualWithAccuracy(real number, real expected, real accuracy)
{
    if (number > (expected - accuracy) && number < (expected + accuracy)) return true;
    return false;
}

MATCHER_P2(RealNearForContainer, containerExpected, accuracy, "")
{
    if (arg.size() != containerExpected.size()) {
        std::cout << "The checked container does not have the same size as the expected container.\n" << std::endl;
        return false;
    }

    for (int i = 0; i < arg.size(); i++) {
        if (!isEqualWithAccuracy(arg[i], containerExpected[i], accuracy)) {
            std::cout << "First mismatching element at index " << i << ": The actual element " << std::to_string(arg[i])
                      << " is not near the expected element " << std::to_string(containerExpected[i]) << ".\n";
            return false;
        }
    }
    return true;
}

MATCHER_P2(NearForContainer, containerExpected, matcher, "")
{
    if (arg.size() != containerExpected.size()) {
        std::cout << "The checked container does not have the same size as the expected container.\n" << std::endl;
        return false;
    }

    for (int i = 0; i < arg.size(); i++) {
        testing::ExplainMatchResult(matcher, arg[i], result_listener);
    }
    return true;
}


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

class CoordinateTransformationTest : public testing::Test
{
protected:
    real datumX = 1.;
    real datumY = 2.;
    real datumZ = 3.;
    std::array<real, 3> data = { 1., 2., 3. };
    std::array<real, 3> dataBeforeTransformation = data;
};

TEST_F(CoordinateTransformationTest, rotatingToGlobal_zeroRotationAndTranslation)
{
    const std::array<real, 3> angles = { 0.0, 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, testing::Eq(dataBeforeTransformation));
}


TEST_F(CoordinateTransformationTest, rotatingToGlobal_twoPiRotation)
{
    const std::array<real, 3> angles = { c2Pi, 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piRotation)
{
    const std::array<real, 3> angles = { cPi, 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, -datumY, -datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piHalfRotation)
{
    const std::array<real, 3> angles = { real(0.5 * cPi), 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, -datumZ, datumY };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_minusPiHalfRotation)
{
    const std::array<real, 3> angles = { real(-0.5 * cPi), 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, datumZ, -datumY };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}
