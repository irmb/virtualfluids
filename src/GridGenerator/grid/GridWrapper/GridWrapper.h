#ifndef GridKernel_H
#define GridKernel_H

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <string>
#include <GridGenerator/grid/Grid.cuh>

struct Vertex;
struct Geometry;

class GridGenerator_EXPORT GridWrapper
{
public:
    virtual ~GridWrapper() {};

    virtual void meshGrid(Geometry &geom) = 0;
    virtual void deleteSolidNodes() = 0;
    virtual void floodFill(const Vertex &vec) = 0;
    virtual void copyDataFromGPU() = 0;

    Grid grid;

private:
	virtual void allocDistribution() = 0;
	virtual void allocField() = 0;
};

#endif
