#ifndef GRID_MOCKS_H
#define GRID_MOCKS_H

#include "GridGenerator/global.h"
#include "GridStrategy/GridStrategyMocks.h"
#include "distributions/Distribution.h"


struct GridDummy
{
private:
    GridDummy(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta) {}

public:
    static SPtr<GridDummy> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d)
    {
        return SPtr<GridDummy>(new GridDummy(startX, startY, startZ, endX, endY, endZ, delta));
    }

    real getDelta() const
    {
        return delta;
    }

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX = 0.0, endY = 0.0, endZ = 0.0;

    real delta;
};

#endif
