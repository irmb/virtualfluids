#ifndef Geometry_h
#define Geometry_h

#include "GridGenerator/global.h"


#include <stdio.h>
#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <memory>

#include "../Triangle/Triangle.h"
#include "../BoundingBox/BoundingBox.h"

#include "../Object.h"

class GeometryMemento;

enum class DiscretizationMethod { RAYCASTING, POINT_IN_OBJECT, POINT_UNDER_TRIANGLE };


struct TriangularMesh : public Object
{
public:

    VF_PUBLIC static TriangularMesh* make(const std::string& fileName, DiscretizationMethod discretizationMethod);
	VF_PUBLIC TriangularMesh();
    VF_PUBLIC TriangularMesh(const std::string& inputPath);
	VF_PUBLIC TriangularMesh(const std::string& inputPath, const BoundingBox &box);
	HOSTDEVICE VF_PUBLIC ~TriangularMesh();

	VF_PUBLIC void setTriangles(std::vector<Triangle> triangles);
	VF_PUBLIC void setMinMax(BoundingBox minmax);

	std::vector<Triangle> triangleVec;
	Triangle *triangles;
	int size;
	BoundingBox minmax;

    HOST VF_PUBLIC bool operator==(const TriangularMesh &geometry) const;

    HOST VF_PUBLIC GeometryMemento getState() const;
    HOST VF_PUBLIC void setState(const GeometryMemento &memento);
    VF_PUBLIC void findNeighbors();

    HOST VF_PUBLIC DiscretizationMethod getDiscretizationMethod() const;


private:
	
	void initalizeDataFromTriangles();

    DiscretizationMethod discretizationMethod;

public:
    HOSTDEVICE Object* clone() const override { return new TriangularMesh(); }
    double getX1Centroid() override { throw "Not implemented in TriangularMesh"; }
    double getX1Minimum() override { return minmax.minX; }
    double getX1Maximum() override { return minmax.maxX; }
    double getX2Centroid() override { throw "Not implemented in TriangularMesh"; }
    double getX2Minimum() override { return minmax.minY; }
    double getX2Maximum() override { return minmax.maxY; }
    double getX3Centroid() override { throw "Not implemented in TriangularMesh"; }
    double getX3Minimum() override { return minmax.minZ; }
    double getX3Maximum() override { return minmax.maxZ; }
    void scale(double delta) override { throw "Not implemented in TriangularMesh"; }
    HOSTDEVICE bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset,
        const double& maxOffset) override {
        return false;
    }
};



#endif

