#ifndef GRID_STRATEGY_H
#define GRID_STRATEGY_H

#include <VirtualFluidsDefinitions.h>
#include "core/PointerDefinitions.h"

struct Vertex;
struct Geometry;
class GridImp;

class VF_PUBLIC GridStrategy
{
public:
    virtual ~GridStrategy() {};

    virtual void allocateGridMemory(SPtr<GridImp> grid) = 0;

    virtual void initalNodes(SPtr<GridImp> grid) = 0;
    virtual void mesh(SPtr<GridImp> grid, Geometry &geom) = 0;

    virtual void createGridInterface(SPtr<GridImp> grid, SPtr<GridImp> finerGrid) = 0;

    virtual void deleteSolidNodes(SPtr<GridImp> grid) = 0;

    virtual void freeMemory(SPtr<GridImp> grid) = 0;

};

#endif
