#ifndef GRID_CPU_STRATEGY_H
#define GRID_CPU_STRATEGY_H

#include "GridGenerator/global.h"

#include "../GridStrategy.h"

#include "core/PointerDefinitions.h"


class GridImp;
class TriangularMesh;

class VF_PUBLIC GridCpuStrategy : public GridStrategy
{
public:
    virtual ~GridCpuStrategy() {};

    void allocateGridMemory(SPtr<GridImp> grid) override;

    void initalNodesToOutOfGrid(SPtr<GridImp> grid) override;
    void findInnerNodes(SPtr<GridImp> grid) override;
    void findInnerNodes(SPtr<GridImp> grid, TriangularMesh* triangularMesh) override;
    void findStopperNodes(SPtr<GridImp> grid) override;

    void mesh(SPtr<GridImp> grid, TriangularMesh &geom) override;

    void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid) override;

    void freeMemory(SPtr<GridImp> grid) override;


    void deleteSolidNodes(SPtr<GridImp> grid) override;

    virtual void copyDataFromGPU() {};

protected:
    static void findForNeighborsNewIndices(SPtr<GridImp> grid);
    static void findForGridInterfaceNewIndices(SPtr<GridImp> grid, SPtr<GridImp> fineGrid);

    void pointUnderTriangleMethod(SPtr<GridImp> grid, TriangularMesh* triangularMesh);
    void rayCastingMethod(SPtr<GridImp> grid, TriangularMesh* triangularMesh);
    void pointInObjectMethod(SPtr<GridImp> grid, TriangularMesh* triangularMesh);
    void removeOddBoundaryCellNodes(SPtr<GridImp> grid);
public:
    void allocateFieldMemory(Field* field) override;
    void freeFieldMemory(Field* field) override;
    void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override;

};

#endif