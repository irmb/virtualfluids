#ifndef BoundingBox_h
#define BoundingBox_h

#include "GridGenerator/global.h"


#include <vector>
#include "cuda_runtime.h"

struct Vertex;
struct Triangle;

class BoundingBoxMemento;

template <typename T>
class VF_PUBLIC BoundingBox
{
public:
	T minX;
	T maxX;
	T minY;
	T maxY;
	T minZ;
	T maxZ;

	BoundingBox(T minX, T maxX, T minY, T maxY, T minZ, T maxZ);
    HOSTDEVICE BoundingBox();
	BoundingBox(const BoundingBox<int> &box);
	BoundingBox(const BoundingBox<real> &t);

public:
	HOST static BoundingBox<real> makeExactBox(const Triangle &t);
	HOSTDEVICE static BoundingBox<int> makeNodeBox(const Triangle &t);
    HOST static BoundingBox<T> makeInvalidMinMaxBox();

    void setMinMax(const Triangle& t);
	void print() const;

	bool isInside(const Triangle &t) const;
	bool intersect(const Triangle &t) const;

	std::vector<std::vector<Vertex> > getIntersectionPoints(const BoundingBox<real> &b) const;
	bool intersect(const BoundingBox<T> &box) const;


    HOST BoundingBoxMemento getState() const;
    HOST void setState(const BoundingBoxMemento &memento);

    HOST bool operator==(const BoundingBox<real> &box) const;
    

private:
    HOSTDEVICE static void calculateMinMaxOnNodes(int &minNode, int &maxNode, const real &minExact, const real &maxExact);

	bool isInside(const Vertex &v) const;
	void getPoints(Vertex v[8]) const;

};

#endif
