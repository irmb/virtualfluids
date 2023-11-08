#include "GeometryUtils.h"
#include "basics/constants/NumericConstants.h"
#include "tests/testUtilities.h"

TEST(GeometryUtilsTest, translate2D)
{
    real newPositionX;
    real newPositionY;

    translate2D(0., 0., newPositionX, newPositionY, 0., 0.);
    EXPECT_THAT(newPositionX, RealEq(0.));
    EXPECT_THAT(newPositionY, RealEq(0.));

    translate2D(0.5, 0.5, newPositionX, newPositionY, 1., 1.);
    EXPECT_THAT(newPositionX, RealEq(1.5));
    EXPECT_THAT(newPositionY, RealEq(1.5));

    translate2D(0.5, 0.5, newPositionX, newPositionY, -1., -1.);
    EXPECT_THAT(newPositionX, RealEq(-0.5));
    EXPECT_THAT(newPositionY, RealEq(-0.5));
}

TEST(GeometryUtilsTest, inverseTranslate2D)
{
    real newPositionX;
    real newPositionY;

    invTranslate2D(0., 0., newPositionX, newPositionY, 0., 0.);
    EXPECT_THAT(newPositionX, RealEq(0.));
    EXPECT_THAT(newPositionY, RealEq(0.));

    invTranslate2D(0.5, 0.5, newPositionX, newPositionY, 1., 1.);
    EXPECT_THAT(newPositionX, RealEq(-0.5));
    EXPECT_THAT(newPositionY, RealEq(-0.5));

    invTranslate2D(0.5, 0.5, newPositionX, newPositionY, -1., -1.);
    EXPECT_THAT(newPositionX, RealEq(1.5));
    EXPECT_THAT(newPositionY, RealEq(1.5));
}

TEST(GeometryUtilsTest, rotate2dAround0)
{
    auto posX = 2.0;
    auto posY = 0.0;
    real newPosX;
    real newPosY;

    auto angle = 0.0;
    rotate2D(angle, posX, posY, newPosX, newPosY);
    EXPECT_THAT(newPosX, RealNear(2.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(0.0, 10e-5));

    angle = 0.5 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY);
    EXPECT_THAT(newPosX, RealNear(0.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(2.0, 10e-5));

    angle = 1.0 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY);
    EXPECT_THAT(newPosX, RealNear(-2.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(0.0, 10e-5));

    angle = 1.5 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY);
    EXPECT_THAT(newPosX, RealNear(0.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-2.0, 10e-5));

    angle = 2.0 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY);
    EXPECT_THAT(newPosX, RealNear(2.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(0.0, 10e-5));

    angle = -0.5 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY);
    EXPECT_THAT(newPosX, RealNear(0.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-2.0, 10e-5));
}

TEST(GeometryUtilsTest, rotate2dWithOrigin)
{
    auto posX = 3.0;
    auto posY = -1.0;
    auto originX = 1.0;
    auto originY = -1.0;
    real newPosX;
    real newPosY;

    auto angle = 0.0;
    rotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));

    angle = 0.5 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));

    angle = 1.0 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(-1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));

    angle = 1.5 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-3.0, 10e-5));

    angle = 2.0 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));

    angle = -0.5 * vf::basics::constant::cPi;
    rotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-3.0, 10e-5));
}

TEST(GeometryUtilsTest, inverseRotate2DWithOrigin)
{
    auto posX = 3.0;
    auto posY = -1.0;
    auto originX = 1.0;
    auto originY = -1.0;
    real newPosX;
    real newPosY;

    auto angle = 0.0;
    invRotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));

    angle = 0.5 * vf::basics::constant::cPi;
    invRotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-3.0, 10e-5));

    angle = 1.0 * vf::basics::constant::cPi;
    invRotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(-1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));

    angle = 1.5 * vf::basics::constant::cPi;
    invRotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));

    angle = 2.0 * vf::basics::constant::cPi;
    invRotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));

    angle = -1.5 * vf::basics::constant::cPi;
    invRotate2D(angle, posX, posY, newPosX, newPosY, originX, originY);
    EXPECT_THAT(newPosX, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-3.0, 10e-5));
}

TEST(GeometryUtilsTest, rotateAboutX3dAround0)
{
    auto posX = 0.5;
    auto posY = 2.0;
    auto posZ = 0.0;
    real newPosX;
    real newPosY;
    real newPosZ;

    auto angle = 0.0;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(2.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(0.0, 10e-5));

    angle = 0.5 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(0.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(2.0, 10e-5));

    angle = 1.0 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-2.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(0.0, 10e-5));

    angle = 1.5 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(0.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-2.0, 10e-5));

    angle = 2.0 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(2.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(0.0, 10e-5));

    angle = -0.5 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(0.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-2.0, 10e-5));
}

TEST(GeometryUtilsTest, rotateAboutX3dWithOrigin)
{
    auto posX = 0.5;
    auto posY = 3.0;
    auto posZ = -1.0;
    auto originX = -0.75;
    auto originY = 1.0;
    auto originZ = -1.0;
    real newPosX;
    real newPosY;
    real newPosZ;

    auto angle = 0.0;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-1.0, 10e-5));

    angle = 0.5 * vf::basics::constant::cPi;    
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(1.0, 10e-5));

    angle = 1.0 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-1.0, 10e-5));

    angle = 1.5 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-3.0, 10e-5));

    angle = 2.0 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-1.0, 10e-5));

    angle = -0.5 * vf::basics::constant::cPi;
    rotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-3.0, 10e-5));
}

TEST(GeometryUtilsTest, inverseRotateAboutX3dWithOrigin)
{
    auto posX = 0.5;
    auto posY = 3.0;
    auto posZ = -1.0;
    auto originX = -0.75;
    auto originY = 1.0;
    auto originZ = -1.0;
    real newPosX;
    real newPosY;
    real newPosZ;

    auto angle = 0.0;
    invRotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-1.0, 10e-5));

    angle = 0.5 * vf::basics::constant::cPi;
    invRotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-3.0, 10e-5));

    angle = 1.0 * vf::basics::constant::cPi;
    invRotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(-1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-1.0, 10e-5));

    angle = 1.5 * vf::basics::constant::cPi;    
    invRotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(1.0, 10e-5));

    angle = 2.0 * vf::basics::constant::cPi;
    invRotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(3.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(-1.0, 10e-5));

    angle = -0.5 * vf::basics::constant::cPi;
    invRotateAboutX3D(angle, posX, posY, posZ, newPosX, newPosY, newPosZ, originX, originY, originZ);
    EXPECT_THAT(newPosX, RealNear(0.5, 10e-5));
    EXPECT_THAT(newPosY, RealNear(1.0, 10e-5));
    EXPECT_THAT(newPosZ, RealNear(1.0, 10e-5));
}