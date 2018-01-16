#ifndef GRID_GPU_STRATEGY_H
#define GRID_GPU_STRATEGY_H

#include <VirtualFluidsDefinitions.h>

#include "../GridStrategy.h"

template <class T>
class BoundingBox;
struct Geometry;

class VF_PUBLIC GridGpuStrategy : public GridStrategy
{
public:
    virtual ~GridGpuStrategy() {};

    void allocateGridMemory(SPtr<Grid> grid) override;

    void initalNodes(SPtr<Grid> grid) override;
    void mesh(SPtr<Grid> grid, Geometry &geom) override;

    void freeMemory(SPtr<Grid> grid) override;


    void deleteSolidNodes(SPtr<Grid> grid) override;


    virtual void copyDataFromGPU(SPtr<Grid> grid);

	//void markNodesToDeleteOutsideOfGeometry();

private:
    void allocField(SPtr<Grid> grid);
    void allocDistribution(SPtr<Grid> grid);
    void allocNeighborsIndices(SPtr<Grid> grid);
    void allocMatrixIndicesOnGPU(SPtr<Grid> grid);

    void allocAndCopyTrianglesToGPU(Geometry &geom);
    void freeTrianglesFromGPU(const Geometry &geom);


    void allocAndCopyMatrixIndicesToGPU(SPtr<Grid> grid);

    void allocAndCopyFieldToGPU(SPtr<Grid> grid);

    void copyAndFreeFieldFromGPU(SPtr<Grid> grid);
    void copyAndFreeDistributiondFromGPU(SPtr<Grid> grid);

    void copyAndFreeNeighborsToCPU(SPtr<Grid> grid);
    void copyAndFreeMatrixIndicesFromGPU(SPtr<Grid> grid, int size);
};

#endif
