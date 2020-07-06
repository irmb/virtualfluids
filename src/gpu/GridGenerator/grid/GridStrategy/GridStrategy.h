#ifndef GRID_STRATEGY_H
#define GRID_STRATEGY_H

#include "Core/LbmOrGks.h"

#include "global.h"

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

	virtual void allocateQs(SPtr<GridImp> grid) = 0;

    virtual void initalNodesToOutOfGrid(SPtr<GridImp> grid) = 0;
    virtual void fixOddCells(SPtr<GridImp> grid) = 0;
    virtual void findInnerNodes(SPtr<GridImp> grid) = 0;
    virtual void addOverlap(SPtr<GridImp> grid) = 0;

    virtual void fixRefinementIntoWall(SPtr<GridImp> grid) = 0;
    virtual void findStopperNodes(SPtr<GridImp> grid) = 0;
	virtual void findBoundarySolidNodes(SPtr<GridImp> grid) = 0; 
	virtual void findEndOfGridStopperNodes(SPtr<GridImp> grid) = 0;
	virtual void findSolidStopperNodes(SPtr<GridImp> grid) = 0;

    virtual void mesh(SPtr<GridImp> grid, TriangularMesh &geom) = 0;

    virtual uint closeNeedleCells(SPtr<GridImp> grid) = 0;
    virtual uint closeNeedleCellsThinWall(SPtr<GridImp> grid) = 0;

    virtual void findQs(SPtr<GridImp> grid, TriangularMesh &geom) = 0;

    virtual void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid, LbmOrGks lbmOrGks) = 0;

    virtual void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) = 0;


    virtual void freeMemory(SPtr<GridImp> grid) = 0;

};

#endif
