#ifndef GRID_STRATEGY_H
#define GRID_STRATEGY_H

#include <VirtualFluidsDefinitions.h>
#include "core/PointerDefinitions.h"
#include "grid/Field.h"

struct Vertex;
struct Geometry;
class GridImp;

class VF_PUBLIC GridStrategy
{
public:
    virtual ~GridStrategy() {}
    virtual void allocateFieldMemory(Field* field) = 0;
    virtual void freeFieldMemory(Field* field) = 0;

    virtual void allocateGridMemory(SPtr<GridImp> grid) = 0;

    virtual void initalNodes(SPtr<GridImp> grid) = 0;
    virtual void mesh(SPtr<GridImp> grid, Geometry &geom) = 0;

    virtual void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> finerGrid) = 0;

    virtual void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) = 0;

    virtual void deleteSolidNodes(SPtr<GridImp> grid) = 0;

    virtual void freeMemory(SPtr<GridImp> grid) = 0;

};

#endif
