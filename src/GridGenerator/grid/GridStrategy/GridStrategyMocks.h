#ifndef GRID_STRATEGYMOCKS_H
#define GRID_STRATEGYMOCKS_H

#include <VirtualFluidsDefinitions.h>
#include "core/PointerDefinitions.h"

#include "GridStrategy.h"

struct Geometry;
class GridImp;

class VF_PUBLIC GridStrategyDummy : public GridStrategy
{
public:
    virtual ~GridStrategyDummy() {}

    virtual void allocateGridMemory(SPtr<GridImp> grid) override {}

    virtual void initalNodes(SPtr<GridImp> grid) override {}
    virtual void mesh(SPtr<GridImp> grid, Geometry &geom) override {}

    virtual void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> finerGrid) override {}

    virtual void deleteSolidNodes(SPtr<GridImp> grid) override {}

    virtual void freeMemory(SPtr<GridImp> grid) override {}
    virtual void allocateFieldMemory(Field* field) override {}
    virtual void freeFieldMemory(Field* field) override {}

    virtual void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override {}
};

#endif
