#ifndef GridKernelCPU_H
#define GridKernelCPU_H

#include "GridGenerator_EXPORT.h"
#include "GridGenerator/global.h"
#include "../GridWrapper.h"

#include <string>


struct Triangle;
struct Vertex;

template <class T>
class BoundingBox;

class GridGenerator_EXPORT GridWrapperCPU : public GridWrapper
{
public:
    GridWrapperCPU(){};
    GridWrapperCPU(BoundingBox<int> &channel, std::string d3Qxx);
	virtual ~GridWrapperCPU();

    virtual void meshGrid(Geometry &geom);
    virtual void deleteSolidNodes();
    virtual void floodFill(const Vertex &vec);
    virtual void copyDataFromGPU() {};

private:
	doubflo initalUniformGrid3d();
	void initialGridNodes();
	virtual void allocDistribution();
	virtual void allocField();
	void runMeshing(const Geometry &geom);

    void findInvalidNodes();
    void findNeighborIndices();
};

#endif
