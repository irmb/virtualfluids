#ifndef GridKernelGPU_H
#define GridKernelGPU_H

#include <GridGenerator_EXPORT.h>

#include <string>
#include "../GridWrapper.h"

template <class T>
class BoundingBox;
struct Geometry;

class GridGenerator_EXPORT GridWrapperGPU : public GridWrapper
{
public:
    GridWrapperGPU(){};
    GridWrapperGPU(BoundingBox<int> &box, std::string d3Qxx);
    virtual ~GridWrapperGPU();

    virtual void meshGrid(Geometry &geom);
    virtual void deleteSolidNodes();
    virtual void floodFill(const Vertex &vec);
    virtual void copyDataFromGPU();

	void markNodesToDeleteOutsideOfGeometry();

private:
	void initalField(BoundingBox<int> &box, std::string direction);
	void allocField();
	void allocDistribution();
    void allocNeighborsIndices();
    void allocAndCopyMatrixIndicesToGPU();
    void allocMatrixIndicesOnGPU();

    void allocAndCopyTrianglesToGPU(Geometry &geom);
    void allocAndCopyFieldToGPU();

    void copyAndFreeFieldFromGPU();
    void copyAndFreeDistributiondFromGPU();
    void freeTrianglesFromGPU(const Geometry &geom);

    void copyAndFreeNeighborsToCPU();
    void copyAndFreeMatrixIndicesFromGPU(int size);
};

#endif
