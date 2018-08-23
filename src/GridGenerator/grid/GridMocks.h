#ifndef GRID_MOCKS_H
#define GRID_MOCKS_H

#include "GridGenerator/global.h"
#include "GridStrategy/GridStrategyMocks.h"
#include "distributions/Distribution.h"
#include <geometries/Vertex/Vertex.h>

#include "Grid.h"
#include "geometries/Object.h"

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

    virtual const Object* getObject() const { return nullptr; }

    virtual uint getSparseSize() const override { return 0; }
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
    virtual int getSparseIndex(uint matrixIndex) const override { return 0; }
    virtual char getFieldEntry(uint matrixIndex) const override { return 0; }
    virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const override {}
    virtual uint* getCF_coarse() const override { return 0; }
    virtual uint* getCF_fine() const override { return 0; }
    virtual uint* getCF_offset() const override { return 0; }
    virtual uint* getFC_coarse() const override { return 0; }
    virtual uint* getFC_fine() const override { return 0; }
    virtual uint* getFC_offset() const override { return 0; }
    virtual real* getDistribution() const override { return nullptr; }
    virtual int* getDirection() const override { return nullptr; }
    virtual int getStartDirection() const override { return 0; }
    virtual int getEndDirection() const override { return 0; }
    virtual void getNodeValues(real* xCoords, real* yCoords, real* zCoords, unsigned* neighborX, unsigned* neighborY,
        unsigned* neighborZ, unsigned* geo) const override {}
    virtual SPtr<GridStrategy> getGridStrategy() const override { return nullptr; }
    virtual void transIndexToCoords(uint index, real& x, real& y, real& z) const override {}
    virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override {}
    virtual void freeMemory() override {}

    virtual void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) override {}
    virtual void mesh(TriangularMesh& geometry) override {}
    virtual uint transCoordToIndex(const real& x, const real& y, const real& z) const override { return 0; }
    virtual int* getNeighborsX() const override { return nullptr; }
    virtual int* getNeighborsY() const override { return nullptr; }
    virtual int* getNeighborsZ() const override { return nullptr; }
    virtual void inital(const SPtr<Grid> fineGrid, uint numberOfLayers) override {};
    virtual void setOddStart( bool xOddStart, bool yOddStart, bool zOddStart ) {};
    virtual bool nodeInCellIs(Cell& cell, char type) const override { return false; }
    virtual void findSparseIndices(SPtr<Grid> fineGrid) override {}
    virtual Vertex getMinimumOnNode(Vertex exact) const override { return Vertex(0, 0, 0); }
    virtual Vertex getMaximumOnNode(Vertex exact) const override { return Vertex(0, 0, 0); }
    virtual void mesh(Object* object) override {}
    void setPeriodicityX(bool periodicity) override {}
    void setPeriodicityY(bool periodicity) override {}
    void setPeriodicityZ(bool periodicity) override {}
    void findQs(Object* object) override {}
    void setFieldEntry(uint matrixIndex, char type) override {}
    real getFirstFluidNode(real coords[3], int direction, real startCoord) const override { return 0.0; }
    real getLastFluidNode(real coords[3], int direction, real startCoord) const override { return 0.0; }

	uint getNumberOfSolidBoundaryNodes() const override { return 0; }
	void setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes) override {}

	real getQValue(const uint index, const uint dir) const override { return 0.0; }

    void setInnerRegionFromFinerGrid( bool innerRegionFromFinerGrid ) override {};
};

class GridStub : public GridDummy
{
protected:
    GridStub(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta), level(level) {}

public:
    static SPtr<GridStub> makeShared(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level)
    {
        return SPtr<GridStub>(new GridStub(gridShape, startX, startY, startZ, endX, endY, endZ, delta, level));
    }

    virtual real getDelta() const override { return delta; }

    virtual real getStartX() const override { return startX; }
    virtual real getStartY() const override { return startY; }
    virtual real getStartZ() const override { return startZ; }
    virtual real getEndX() const override { return endX; }
    virtual real getEndY() const override { return endY; }
    virtual real getEndZ() const override { return endZ; }


private:
    real startX, startY, startZ;
    real endX, endY, endZ;

    real delta;

    uint level;
};


class GridSpy : public GridStub
{
protected:
    GridSpy(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) : GridStub(gridShape, startX, startY, startZ, endX, endY, endZ, delta, level) {}

public:
    static SPtr<GridSpy> makeShared(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level)
    {
        return SPtr<GridSpy>(new GridSpy(gridShape, startX, startY, startZ, endX, endY, endZ, delta, level));
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

    virtual void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) override
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
