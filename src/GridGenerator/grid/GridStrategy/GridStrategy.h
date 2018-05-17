#ifndef GRID_STRATEGY_H
#define GRID_STRATEGY_H

#include <VirtualFluidsDefinitions.h>
#include "core/PointerDefinitions.h"
#include "grid/Field.h"

struct Vertex;
class TriangularMesh;
class GridImp;

class VF_PUBLIC GridStrategy
{
public:
    virtual ~GridStrategy() {}
    virtual void allocateFieldMemory(Field* field) = 0;
    virtual void freeFieldMemory(Field* field) = 0;

    virtual void allocateGridMemory(SPtr<GridImp> grid) = 0;

    virtual void initalNodesToOutOfGrid(SPtr<GridImp> grid) = 0;

    virtual void findInnerNodes(SPtr<GridImp> grid) = 0;
    virtual void findStopperNodes(SPtr<GridImp> grid) = 0;

    virtual void mesh(SPtr<GridImp> grid, TriangularMesh &geom) = 0;

    virtual void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> finerGrid) = 0;

    virtual void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) = 0;


    virtual void freeMemory(SPtr<GridImp> grid) = 0;

};

#endif
