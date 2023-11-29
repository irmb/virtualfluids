#ifndef GRIDDIMENSIONS_H
#define GRIDDIMENSIONS_H

#include "DataTypes.h"

struct GridDimensions {
public:
    real minX;
    real maxX;
    real minY;
    real maxY;
    real minZ;
    real maxZ;
    real delta;

    GridDimensions() = default;
    GridDimensions(real minX, real maxX, real minY, real maxY, real minZ, real maxZ, real delta)
        : minX(minX), maxX(maxX), minY(minY), maxY(maxY), minZ(minZ), maxZ(maxZ), delta(delta){};
};

#endif