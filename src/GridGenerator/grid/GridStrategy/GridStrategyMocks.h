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
    virtual ~GridStrategyDummy() {};

    virtual void allocateGridMemory(SPtr<GridImp> grid) override {};

    virtual void initalNodes(SPtr<GridImp> grid) override {};
    virtual void mesh(SPtr<GridImp> grid, Geometry &geom) override {};

    virtual void createGridInterface(SPtr<GridImp> grid, SPtr<GridImp> finerGrid) override {};

    virtual void deleteSolidNodes(SPtr<GridImp> grid) override {};

    virtual void freeMemory(SPtr<GridImp> grid) override {};

};

#endif
