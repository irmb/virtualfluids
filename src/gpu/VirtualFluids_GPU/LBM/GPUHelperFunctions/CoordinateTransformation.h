#ifndef COORDINATE_TRANSFORMATION
#define COORDINATE_TRANSFORMATION

#include <basics/DataTypes.h>
#include <math.h>
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;

__inline__ __host__ __device__ void rotateDataFromRotatingToGlobal(real &datumX, real &datumY, real &datumZ, real angleX, real angleY,
                                                          real angleZ)
{
    real datumXTemp = datumX;
    real datumYTemp = datumY;
    real datumZTemp = datumZ;
    if (angleX != c0o1) {
        datumYTemp = datumY * cos(angleX) - datumZ * sin(angleX);
        datumZTemp = datumY * sin(angleX) + datumZ * cos(angleX);
    } else if (angleY != c0o1) {
        // rotate in y
        datumXTemp = datumX * cos(angleY) + datumZ * sin(angleY);
        datumZTemp = -datumX * sin(angleY) + datumZ * cos(angleY);
    } else if (angleZ != c0o1) {
        // rotate in z
        datumXTemp = datumX * cos(angleZ) - datumY * sin(angleZ);
        datumYTemp = datumX * sin(angleZ) + datumY * cos(angleZ);
    }
    datumX = datumXTemp;
    datumY = datumYTemp;
    datumZ = datumZTemp;
}

__inline__ __host__ __device__ void transformRotatingToGlobal(real &globalX, real &globalY, real &globalZ, real localX, real localY,
                                                     real localZ, real centerCoordX, real centerCoordY, real centerCoordZ,
                                                     real angleX, real angleY, real angleZ)
{
    globalX = localX;
    globalY = localY;
    globalZ = localZ;

    rotateDataFromRotatingToGlobal(globalX, globalY, globalZ, angleX, angleY, angleZ);

    // translate
    globalX += centerCoordX;
    globalY += centerCoordY;
    globalZ += centerCoordZ;
}


__inline__ __host__ __device__ void rotateDataFromGlobalToRotating(real &datumX, real &datumY, real &datumZ, real angleX,
                                                              real angleY, real angleZ)
{
    real datumXTemp = datumX;
    real datumYTemp = datumY;
    real datumZTemp = datumZ;
    if (angleX != c0o1) {
        datumYTemp = datumY * cos(angleX) + datumZ * sin(angleX);
        datumZTemp = -datumY * sin(angleX) + datumZ * cos(angleX);
    } else if (angleY != c0o1) {
        // rotate in y
        datumXTemp = datumX * cos(angleY) - datumZ * sin(angleY);
        datumZTemp = datumX * sin(angleY) + datumZ * cos(angleY);
    } else if (angleZ != c0o1) {
        // rotate in z
        datumXTemp = datumX * cos(angleZ) + datumY * sin(angleZ);
        datumYTemp = -datumX * sin(angleZ) + datumY * cos(angleZ);
    }
    datumX = datumXTemp;
    datumY = datumYTemp;
    datumZ = datumZTemp;
}

__inline__ __host__ __device__ void transformGlobalToRotating(real &rotatingX, real &rotatingY, real &rotatingZ, real globalX,
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

    rotateDataFromGlobalToRotating(rotatingX, rotatingY, rotatingZ, angleX, angleY, angleZ);
}

#endif