#ifndef COORDINATE_TRANSFORMATION
#define COORDINATE_TRANSFORMATION

#include <basics/DataTypes.h>
#include <math.h>
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;

__inline__ __device__ void rotateVelocityFromRotatingToGlobal(real &velocityX, real &velocityY, real &velocityZ, real angleX,
                                                              real angleY, real angleZ)
{
    real velocityXTemp = velocityX;
    real velocityYTemp = velocityY;
    real velocityZTemp = velocityZ;
    if (angleX != c0o1) {
        velocityYTemp = velocityY * cos(angleX) - velocityZ * sin(angleX);
        velocityZTemp = velocityY * sin(angleX) + velocityZ * cos(angleX);
    } else if (angleY != c0o1) {
        // rotate in y
        velocityXTemp = velocityX * cos(angleY) + velocityZ * sin(angleY);
        velocityZTemp = -velocityX * sin(angleY) + velocityZ * cos(angleY);
    } else if (angleZ != c0o1) {
        // rotate in z
        velocityXTemp = velocityX * cos(angleZ) - velocityY * sin(angleZ);
        velocityYTemp = velocityX * sin(angleZ) + velocityY * cos(angleZ);
    }
    velocityX = velocityXTemp;
    velocityY = velocityYTemp;
    velocityZ = velocityZTemp;
}

__inline__ __device__ void transformRotatingToGlobal(real &globalX, real &globalY, real &globalZ, real localX, real localY,
                                                     real localZ, real centerCoordX, real centerCoordY, real centerCoordZ,
                                                     real angleX, real angleY, real angleZ)
{
    globalX = localX;
    globalY = localY;
    globalZ = localZ;

    // rotate
    if (angleX != 0) {
        // rotate in x
        globalY = localY * cos(angleX) - localZ * sin(angleX);
        globalZ = localY * sin(angleX) + localZ * cos(angleX);
    } else if (angleY != 0) {
        // rotate in y
        globalX = localX * cos(angleY) + localZ * sin(angleY);
        globalZ = -localX * sin(angleY) + localZ * cos(angleY);
    } else if (angleZ != 0) {
        // rotate in z
        globalX = localX * cos(angleZ) - localY * sin(angleZ);
        globalY = localX * sin(angleZ) + localY * cos(angleZ);
    }

    // translate
    globalX += centerCoordX;
    globalY += centerCoordY;
    globalZ += centerCoordZ;
}


__inline__ __device__ void rotateVelocityFromGlobalToRotating(real &velocityX, real &velocityY, real &velocityZ, real angleX,
                                                              real angleY, real angleZ)
{
    real velocityXTemp = velocityX;
    real velocityYTemp = velocityY;
    real velocityZTemp = velocityZ;
    if (angleX != c0o1) {
        velocityYTemp = velocityY * cos(angleX) + velocityZ * sin(angleX);
        velocityZTemp = -velocityY * sin(angleX) + velocityZ * cos(angleX);
    } else if (angleY != c0o1) {
        // rotate in y
        velocityXTemp = velocityX * cos(angleY) - velocityZ * sin(angleY);
        velocityZTemp = velocityX * sin(angleY) + velocityZ * cos(angleY);
    } else if (angleZ != c0o1) {
        // rotate in z
        velocityXTemp = velocityX * cos(angleZ) + velocityY * sin(angleZ);
        velocityYTemp = -velocityX * sin(angleZ) + velocityY * cos(angleZ);
    }
    velocityX = velocityXTemp;
    velocityY = velocityYTemp;
    velocityZ = velocityZTemp;
}

__inline__ __device__ void transformGlobalToRotating(real &rotatingX, real &rotatingY, real &rotatingZ, real globalX,
                                                              real globalY, real globalZ, real centerCoordX, real centerCoordY,
                                                              real centerCoordZ, real angleX, real angleY, real angleZ)
{

    // translate
    globalX -= centerCoordX;
    globalY -= centerCoordY;
    globalZ -= centerCoordZ;

    rotatingX = globalX;
    rotatingY = globalY;
    rotatingZ = globalZ;

    // rotate
    if (angleX != 0) {
        // rotate in x
        rotatingY = globalY * cos(angleX) + globalZ * sin(angleX);
        rotatingZ = -globalY * sin(angleX) + globalZ * cos(angleX);
    } else if (angleY != 0) {
        // rotate in y
        rotatingX = globalX * cos(angleY) - globalZ * sin(angleY);
        rotatingZ = globalX * sin(angleY) + globalZ * cos(angleY);
    } else if (angleZ != 0) {
        // rotate in z
        rotatingX = globalX * cos(angleZ) + globalY * sin(angleZ);
        rotatingY = -globalX * sin(angleZ) + globalY * cos(angleZ);
    }
}

#endif