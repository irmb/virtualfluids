#ifndef COORDINATE_TRANSFORMATION
#define COORDINATE_TRANSFORMATION

#include <basics/DataTypes.h>
#include <math.h>

__inline__ __device__ void transformRotatingToGlobal(real &globalX, real &globalY, real &globalZ, real localX,
                                                              real localY, real localZ, real centerCoordX, real centerCoordY,
                                                              real centerCoordZ, real angleX, real angleY, real angleZ)
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

#endif