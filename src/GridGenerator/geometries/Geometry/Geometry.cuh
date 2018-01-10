#ifndef Geometry_h
#define Geometry_h

#include "GridGenerator/global.h"


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
	VF_PUBLIC Geometry();
	VF_PUBLIC Geometry(const Geometry& geo);
	VF_PUBLIC Geometry(const std::string& inputPath, const BoundingBox<int> &box, const Transformator *trafo);
	VF_PUBLIC ~Geometry();
	VF_PUBLIC Transformator* getTransformator();

	VF_PUBLIC void setTriangles(std::vector<Triangle> triangles);
	VF_PUBLIC void setMinMax(BoundingBox<real> minmax);
	VF_PUBLIC void transformChannelGeometry(const real resolution);

	std::vector<Triangle> triangleVec;
	Triangle *triangles;
	int size;
	BoundingBox<real> minmax;

    HOST VF_PUBLIC bool operator==(const Geometry &geometry) const;

    HOST VF_PUBLIC GeometryMemento getState() const;
    HOST VF_PUBLIC void setState(const GeometryMemento &memento);



private:
	
	Transformator* transformator;
	void findNeighbors();
	void initalizeDataFromTriangles();

};



#endif

