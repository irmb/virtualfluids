#ifndef BoundingBox_h
#define BoundingBox_h

#include <vector>
#include <cuda_runtime.h>

#include "global.h"

struct Vertex;
struct Triangle;


class VIRTUALFLUIDS_GPU_EXPORT BoundingBox
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
    CUDA_HOST static BoundingBox makeInvalidMinMaxBox();

    void setMinMax(const Triangle& t);
	void print() const;

	bool isInside(const Triangle &t) const;
	bool isInside(const real x, const real y, const real z) const;
	bool intersect(const Triangle &t) const;

	std::vector<std::vector<Vertex> > getIntersectionPoints(const BoundingBox &b) const;
	bool intersect(const BoundingBox &box) const;

    CUDA_HOST bool operator==(const BoundingBox &box) const;
    
    void extend(real delta);

private:
	bool isInside(const Vertex &v) const;
	void getPoints(Vertex v[8]) const;

};

#endif

