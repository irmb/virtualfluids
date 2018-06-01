#ifndef GRID_STRATEGYMOCKS_H
#define GRID_STRATEGYMOCKS_H

#include <VirtualFluidsDefinitions.h>
#include "core/PointerDefinitions.h"

#include "GridStrategy.h"

class TriangularMesh;
class GridImp;

class VF_PUBLIC GridStrategyDummy : public GridStrategy
{
public:
    virtual ~GridStrategyDummy() {}

    virtual void allocateGridMemory(SPtr<GridImp> grid) override {}

    virtual void initalNodesToOutOfGrid(SPtr<GridImp> grid) override {}
    virtual void findInnerNodes(SPtr<GridImp> grid) override {}
    virtual void findStopperNodes(SPtr<GridImp> grid) override {}
	void findEndOfGridStopperNodes(SPtr<GridImp> grid) override {}
	void findSolidStopperNodes(SPtr<GridImp> grid) override {}

    virtual void mesh(SPtr<GridImp> grid, TriangularMesh &geom) override {}

    virtual void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> finerGrid, LbmOrGks lbmOrGks) override {}


    virtual void freeMemory(SPtr<GridImp> grid) override {}
    virtual void allocateFieldMemory(Field* field) override {}
    virtual void freeFieldMemory(Field* field) override {}

    virtual void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override {}


    void findQs(SPtr<GridImp> grid, TriangularMesh& geom) override {};
};

#endif
