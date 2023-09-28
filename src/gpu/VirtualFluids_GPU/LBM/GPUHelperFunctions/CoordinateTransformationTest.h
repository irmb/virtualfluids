
#include "CoordinateTransformation.h"
#include "tests/testUtilities.h"

using namespace testingVF;

class CoordinateTransformationTest : public testing::Test
{
protected:
    real datumX = 1.;
    real datumY = 2.;
    real datumZ = 3.;
    std::array<real, 3> data = { datumX, datumY, datumZ };
    std::array<real, 3> dataBeforeTransformation = data;
    std::array<real, 3> centerCoordinates = {0.0, 0.0, 0.0};
};

TEST_F(CoordinateTransformationTest, rotatingToGlobal_zeroRotationAndTranslation)
{
    const std::array<real, 3> angles = { 0.0, 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, testing::Eq(dataBeforeTransformation));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, testing::Eq(dataBeforeTransformation));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_twoPiRotationX)
{
    const std::array<real, 3> angles = { c2Pi, 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piRotationX)
{
    const std::array<real, 3> angles = { cPi, 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, -datumY, -datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piHalfRotationX)
{
    const std::array<real, 3> angles = { real(0.5 * cPi), 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, -datumZ, datumY };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_minusPiHalfRotationX)
{
    const std::array<real, 3> angles = { real(-0.5 * cPi), 0.0, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, datumZ, -datumY };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_twoPiRotationY)
{
    const std::array<real, 3> angles = { 0.0, c2Pi, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piRotationY)
{
    const std::array<real, 3> angles = { 0.0, cPi, 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumX, datumY, -datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piHalfRotationY)
{
    const std::array<real, 3> angles = { 0.0, real(0.5 * cPi), 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumZ, datumY, -datumX };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_minusPiHalfRotationY)
{
    const std::array<real, 3> angles = { 0.0, real(-0.5 * cPi), 0.0 };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumZ, datumY, datumX };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_twoPiRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, c2Pi};
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));


    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, cPi };
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumX, -datumY, datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_piHalfRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, real(0.5 * cPi)};
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumY, datumX, datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, rotatingToGlobal_minusPiHalfRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, real(-0.5 * cPi)};
    rotateDataFromRotatingToGlobal(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumY, -datumX, datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformRotatingToGlobal(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}










////////////////////////////////////////////////////////////////////////////////
// global to rotating

TEST_F(CoordinateTransformationTest, globalToRotating_zeroRotationAndTranslation)
{
    const std::array<real, 3> angles = { 0.0, 0.0, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, testing::Eq(dataBeforeTransformation));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, testing::Eq(dataBeforeTransformation));
}

TEST_F(CoordinateTransformationTest, globalToRotating_twoPiRotationX)
{
    const std::array<real, 3> angles = { c2Pi, 0.0, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_piRotationX)
{
    const std::array<real, 3> angles = { cPi, 0.0, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, -datumY, -datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_piHalfRotationX)
{
    const std::array<real, 3> angles = { real(0.5 * cPi), 0.0, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, datumZ, -datumY };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_minusPiHalfRotationX)
{
    const std::array<real, 3> angles = { real(-0.5 * cPi), 0.0, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumX, -datumZ, datumY };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));}

TEST_F(CoordinateTransformationTest, globalToRotating_twoPiRotationY)
{
    const std::array<real, 3> angles = { 0.0, c2Pi, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_piRotationY)
{
    const std::array<real, 3> angles = { 0.0, cPi, 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumX, datumY, -datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));}

TEST_F(CoordinateTransformationTest, globalToRotating_piHalfRotationY)
{
    const std::array<real, 3> angles = { 0.0, real(0.5 * cPi), 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumZ, datumY, datumX };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_minusPiHalfRotationY)
{
    const std::array<real, 3> angles = { 0.0, real(-0.5 * cPi), 0.0 };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumZ, datumY, -datumX };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_twoPiRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, c2Pi};
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));


    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(dataBeforeTransformation, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_piRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, cPi };
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumX, -datumY, datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_piHalfRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, real(0.5 * cPi)};
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { datumY, -datumX, datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}

TEST_F(CoordinateTransformationTest, globalToRotating_minusPiHalfRotationZ)
{
    const std::array<real, 3> angles = { 0.0, 0.0, real(-0.5 * cPi)};
    rotateDataFromGlobalToRotating(data[0], data[1], data[2], angles[0], angles[1], angles[2]);
    const std::array<real, 3> expected = { -datumY, datumX, datumZ };
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));

    transformGlobalToRotating(data[0], data[1], data[2], dataBeforeTransformation[0], dataBeforeTransformation[1],
                              dataBeforeTransformation[2], centerCoordinates[0], centerCoordinates[1], centerCoordinates[2],
                              angles[0], angles[1], angles[2]);
    EXPECT_THAT(data, RealNearForContainer(expected, 1e-6));
}
