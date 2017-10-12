#ifndef Geometry_h
#define Geometry_h

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <stdio.h>
#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <memory>

#include "../Triangle/Triangle.cuh"
#include <GridGenerator/utilities/Transformator/Transformator.h>
#include "../BoundingBox/BoundingBox.cuh"


class GeometryMemento;

struct Geometry
{
public:
	GridGenerator_EXPORT Geometry();
	GridGenerator_EXPORT Geometry(const Geometry& geo);
	GridGenerator_EXPORT Geometry(const std::string& inputPath, const BoundingBox<int> &box, const Transformator *trafo);
	GridGenerator_EXPORT ~Geometry();
	GridGenerator_EXPORT Transformator* getTransformator();

	GridGenerator_EXPORT void setTriangles(std::vector<Triangle> triangles);
	GridGenerator_EXPORT void setMinMax(BoundingBox<doubflo> minmax);
	GridGenerator_EXPORT void transformChannelGeometry(const doubflo resolution);

	std::vector<Triangle> triangleVec;
	Triangle *triangles;
	int size;
	BoundingBox<doubflo> minmax;

    HOST GridGenerator_EXPORT bool operator==(const Geometry &geometry) const;

    HOST GridGenerator_EXPORT GeometryMemento getState() const;
    HOST GridGenerator_EXPORT void setState(const GeometryMemento &memento);



private:
	
	Transformator* transformator;
	void findNeighbors();
	void initalizeDataFromTriangles();

};



#endif

