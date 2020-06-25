#ifndef GRID_CPU_STRATEGY_H
#define GRID_CPU_STRATEGY_H

#include "global.h"

#include "grid/GridStrategy/GridStrategy.h"

class GridImp;
class TriangularMesh;

class VF_PUBLIC GridCpuStrategy : public GridStrategy
{
public:
    virtual ~GridCpuStrategy() {};

    void allocateGridMemory(SPtr<GridImp> grid) override;

	void allocateQs(SPtr<GridImp> grid) override;

    void initalNodesToOutOfGrid(SPtr<GridImp> grid) override;
    void fixOddCells(SPtr<GridImp> grid) override;
    void findInnerNodes(SPtr<GridImp> grid) override;
    void addOverlap(SPtr<GridImp> grid) override;
    void fixRefinementIntoWall(SPtr<GridImp> grid) override;
    void findStopperNodes(SPtr<GridImp> grid) override;
	void findBoundarySolidNodes(SPtr<GridImp> grid)  override;
	void findEndOfGridStopperNodes(SPtr<GridImp> grid) override;
	void findSolidStopperNodes(SPtr<GridImp> grid) override;

    void mesh(SPtr<GridImp> grid, TriangularMesh &geom) override;

    uint closeNeedleCells(SPtr<GridImp> grid) override;
    uint closeNeedleCellsThinWall(SPtr<GridImp> grid) override;

    void findQs(SPtr<GridImp> grid, TriangularMesh &geom) override;

    void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid, LbmOrGks lbmOrGks) override;

    void freeMemory(SPtr<GridImp> grid) override;



    virtual void copyDataFromGPU() {};

protected:
    static void findForNeighborsNewIndices(SPtr<GridImp> grid);
    static void findForGridInterfaceNewIndices(SPtr<GridImp> grid, SPtr<GridImp> fineGrid);
public:
    void allocateFieldMemory(Field* field) override;
    void freeFieldMemory(Field* field) override;
    void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override;

};

#endif