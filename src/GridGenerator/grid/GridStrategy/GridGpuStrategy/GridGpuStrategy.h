#ifndef GRID_GPU_STRATEGY_H
#define GRID_GPU_STRATEGY_H

#include <VirtualFluidsDefinitions.h>
#include "global.h"

#include "../GridStrategy.h"

template <class T>
class BoundingBox;
struct Geometry;

class VF_PUBLIC GridGpuStrategy : public GridStrategy
{
public:
    virtual ~GridGpuStrategy() {};

    void allocateGridMemory(SPtr<GridImp> grid) override;

    void initalNodes(SPtr<GridImp> grid) override;
    void mesh(SPtr<GridImp> grid, Geometry &geom) override;
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

    void allocAndCopyTrianglesToGPU(Geometry &geom);
    void freeTrianglesFromGPU(const Geometry &geom);


    void allocAndCopyMatrixIndicesToGPU(SPtr<GridImp> grid, const uint& size);

    void allocAndCopyFieldToGPU(SPtr<GridImp> grid, const uint& size);

    void copyAndFreeFieldFromGPU(SPtr<GridImp> grid);
    void copyAndFreeDistributiondFromGPU(SPtr<GridImp> grid);

    void copyAndFreeNeighborsToCPU(SPtr<GridImp> grid);
    void copyAndFreeMatrixIndicesFromGPU(SPtr<GridImp> grid, int size);

};

#endif
