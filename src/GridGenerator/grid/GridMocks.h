#ifndef GRID_MOCKS_H
#define GRID_MOCKS_H

#include "GridGenerator/global.h"
#include "GridStrategy/GridStrategyMocks.h"
#include "distributions/Distribution.h"


class GridDummy
{
protected:
    GridDummy() {}

public:
    virtual ~GridDummy() {}

    static SPtr<GridDummy> makeShared()
    {
        return SPtr<GridDummy>(new GridDummy());
    }

    real getDelta() const
    {
        return 0.0;
    }

};

class GridStub : public GridDummy
{
private:
    GridStub(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta) {}

public:
    static SPtr<GridStub> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d)
    {
        return SPtr<GridStub>(new GridStub(startX, startY, startZ, endX, endY, endZ, delta));
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
