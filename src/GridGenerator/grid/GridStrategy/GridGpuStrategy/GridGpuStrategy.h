#ifndef GRID_GPU_STRATEGY_H
#define GRID_GPU_STRATEGY_H

#include <VirtualFluidsDefinitions.h>
#include "global.h"

#include "../GridStrategy.h"

class BoundingBox;
struct TriangularMesh;

class VF_PUBLIC GridGpuStrategy : public GridStrategy
{
public:
    virtual ~GridGpuStrategy() {};

    void allocateGridMemory(SPtr<GridImp> grid) override;

    void initalNodesToOutOfGrid(SPtr<GridImp> grid) override;
    void findInnerNodes(SPtr<GridImp> grid, TriangularMesh* triangularMesh) override;
    void findInnerNodes(SPtr<GridImp> grid) override;
    void findStopperNodes(SPtr<GridImp> grid) override;

    void mesh(SPtr<GridImp> grid, TriangularMesh &geom) override;
    void findGridInterface(SPtr<GridImp> grid, SPtr<GridImp> fineGrid) override;

    void freeMemory(SPtr<GridImp> grid) override;


    void deleteSolidNodes(SPtr<GridImp> grid) override;

    void copyAndFreeGridInterfaceFromGPU(SPtr<GridImp> grid);
    virtual void copyDataFromGPU(SPtr<GridImp> grid);

	//void markNodesToDeleteOutsideOfGeometry();

private:
    void allocField(SPtr<GridImp> grid);
    void allocDistribution(SPtr<GridImp> grid);
    void allocNeighborsIndices(SPtr<GridImp> grid);
    void allocMatrixIndicesOnGPU(SPtr<GridImp> grid);

    void allocAndCopyTrianglesToGPU(TriangularMesh &geom);
    void freeTrianglesFromGPU(const TriangularMesh &geom);


    void allocAndCopyMatrixIndicesToGPU(SPtr<GridImp> grid, const uint& size);

    void allocAndCopyFieldToGPU(Field& field);
    void copyAndFreeFieldFromGPU(Field& field);

    void copyAndFreeDistributiondFromGPU(SPtr<GridImp> grid);

    void copyAndFreeNeighborsToCPU(SPtr<GridImp> grid);
    void copyAndFreeMatrixIndicesFromGPU(SPtr<GridImp> grid, int size);

public:
    void allocateFieldMemory(Field* field) override;
    void freeFieldMemory(Field* field) override;
    void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override;
    
};

#endif
