#ifndef BoundingBox_h
#define BoundingBox_h

#include "GridGenerator/global.h"

#include <vector>
#include "cuda_runtime.h"

struct Vertex;
struct Triangle;


class VF_PUBLIC BoundingBox
{
public:
	real minX;
	real maxX;
	real minY;
	real maxY;
	real minZ;
	real maxZ;

	BoundingBox(real minX, real maxX, real minY, real maxY, real minZ, real maxZ);
    HOSTDEVICE BoundingBox();
	BoundingBox(const BoundingBox &t);

public:
    HOST static BoundingBox makeInvalidMinMaxBox();

    void setMinMax(const Triangle& t);
	void print() const;

	bool isInside(const Triangle &t) const;
	bool intersect(const Triangle &t) const;

	std::vector<std::vector<Vertex> > getIntersectionPoints(const BoundingBox &b) const;
	bool intersect(const BoundingBox &box) const;

    HOST bool operator==(const BoundingBox &box) const;
    

private:
	bool isInside(const Vertex &v) const;
	void getPoints(Vertex v[8]) const;

};

#endif

