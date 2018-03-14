#ifndef GRID_CPU_STRATEGY_H
#define GRID_CPU_STRATEGY_H

#include "GridGenerator/global.h"

#include "../GridStrategy.h"

#include "core/PointerDefinitions.h"


class GridImp;
struct Geometry;

class VF_PUBLIC GridCpuStrategy : public GridStrategy
{
public:
	virtual ~GridCpuStrategy() {};

    void allocateGridMemory(SPtr<GridImp> grid) override;

    void initalNodes(SPtr<GridImp> grid) override;
    void mesh(SPtr<GridImp> grid, Geometry &geom) override;

    void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid) override;

    void freeMemory(SPtr<GridImp> grid) override;


    void deleteSolidNodes(SPtr<GridImp> grid) override;

    virtual void copyDataFromGPU() {};

protected:
    static void findInvalidNodes(SPtr<GridImp> grid);
    static void findForNeighborsNewIndices(SPtr<GridImp> grid);
    static void findForGridInterfaceNewIndices(SPtr<GridImp> grid, SPtr<GridImp> fineGrid);


public:
    void allocateFieldMemory(Field* field) override;
    void freeFieldMemory(Field* field) override;
    void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override;

};

#endif
