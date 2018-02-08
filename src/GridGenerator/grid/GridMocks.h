#ifndef GRID_MOCKS_H
#define GRID_MOCKS_H

#include "GridGenerator/global.h"
#include "GridStrategy/GridStrategyMocks.h"
#include "distributions/Distribution.h"

#include "Grid.h"

class GridDummy : public Grid
{
public:
    virtual ~GridDummy() {}

    static SPtr<GridDummy> makeShared()
    {
        return SPtr<GridDummy>(new GridDummy());
    }

    virtual real getDelta() const
    {
        return 0.0;
    }

    virtual uint getReducedSize() const override { return 0; }
    virtual uint getSize() const override { return 0; }
    virtual real getStartX() const override { return 0.0; }
    virtual real getStartY() const override { return 0.0; }
    virtual real getStartZ() const override { return 0.0; }
    virtual real getEndX() const override { return 0.0; }
    virtual real getEndY() const override { return 0.0; }
    virtual real getEndZ() const override { return 0.0; }
    virtual uint getNumberOfNodesX() const override { return 0; }
    virtual uint getNumberOfNodesY() const override { return 0; }
    virtual uint getNumberOfNodesZ() const override { return 0; }
    virtual uint getNumberOfNodesCF() const override { return 0; }
    virtual uint getNumberOfNodesFC() const override { return 0; }
    virtual int getIndex(uint matrixIndex) const override { return 0; }
    virtual char getFieldEntry(uint matrixIndex) const override { return 0; }
    virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const override {}
    virtual uint* getCF_coarse() const override { return 0; }
    virtual uint* getCF_fine() const override { return 0; }
    virtual uint* getFC_coarse() const override { return 0; }
    virtual uint* getFC_fine() const override { return 0; }
    virtual real* getDistribution() const override { return nullptr; }
    virtual int* getDirection() const override { return nullptr; }
    virtual int getStartDirection() const override { return 0; }
    virtual int getEndDirection() const override { return 0; }
    virtual void setNodeValues(real* xCoords, real* yCoords, real* zCoords, unsigned* neighborX, unsigned* neighborY,
        unsigned* neighborZ, unsigned* geo) const override {}
    virtual SPtr<GridStrategy> getGridStrategy() const override { return nullptr; }
    virtual void transIndexToCoords(int index, real& x, real& y, real& z) const override {}
    virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override {}
    virtual void freeMemory() override {}

    virtual void setStartX(real startX) override {}
    virtual void setStartY(real startY) override {}
    virtual void setStartZ(real startZ) override {}
    virtual void setEndX(real endX) override {}
    virtual void setEndY(real endY) override {}
    virtual void setEndZ(real endZ) override {}
    virtual void removeOverlapNodes(SPtr<Grid> grid) override {}
    virtual void mesh(Geometry& geometry) override {}
    virtual int transCoordToIndex(const Vertex& v) const override { return 0; }
    virtual int transCoordToIndex(const real& x, const real& y, const real& z) const override { return 0; }
    virtual void setFieldEntry(uint index, char entry) override {}
    virtual int* getNeighborsX() const override { return nullptr; }
    virtual int* getNeighborsY() const override { return nullptr; }
    virtual int* getNeighborsZ() const override { return nullptr; }
    virtual void allocateGridMemory() override {}
};

class GridStub : public GridDummy
{
protected:
    GridStub(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta) {}

public:
    static SPtr<GridStub> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d)
    {
        return SPtr<GridStub>(new GridStub(startX, startY, startZ, endX, endY, endZ, delta));
    }

    virtual real getDelta() const override { return delta; }

    virtual real getStartX() const override { return startX; }
    virtual real getStartY() const override { return startY; }
    virtual real getStartZ() const override { return startZ; }
    virtual real getEndX() const override { return endX; }
    virtual real getEndY() const override { return endY; }
    virtual real getEndZ() const override { return endZ; }

    virtual void setStartX(real startX) override { this->startX = startX; }
    virtual void setStartY(real startY) override { this->startY = startY; }
    virtual void setStartZ(real startZ) override { this->startZ = startZ; }
    virtual void setEndX(real endX) override { this->endX = endX; }
    virtual void setEndY(real endY) override { this->endY = endY; }
    virtual void setEndZ(real endZ) override { this->endZ = endZ; }

private:
    real startX, startY, startZ;
    real endX, endY, endZ;

    real delta;
};


class GridSpy : public GridStub
{
private:
    GridSpy(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta) : GridStub(startX, startY, startZ, endX, endY, endZ, delta) {}

public:
    static SPtr<GridSpy> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d)
    {
        return SPtr<GridSpy>(new GridSpy(startX, startY, startZ, endX, endY, endZ, delta));
    }

    bool hasGridInterface() const
    {
        return _hasGridInterface;
    }

    bool hasNoPeriodicityBoundaries() const
    {
        return !(periodicityX && periodicityY && periodicityZ);
    }

    virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override
    {
        this->periodicityX = periodicityX;
        this->periodicityY = periodicityY;
        this->periodicityZ = periodicityZ;
    }

    virtual void removeOverlapNodes(SPtr<Grid> grid) override
    {
        _hasGridInterface = true;
    }

private:
    bool _hasGridInterface = false;
    bool periodicityX;
    bool periodicityY;
    bool periodicityZ;
};

#endif
