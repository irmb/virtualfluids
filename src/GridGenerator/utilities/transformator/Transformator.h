#ifndef Transformator_h
#define Transformator_h

#include <memory>
#include "GridGenerator/global.h"


template<typename T>
class BoundingBox;
struct Triangle;
struct TriangularMesh;
struct Vertex;


class Transformator
{
public:
    static VF_PUBLIC std::shared_ptr<Transformator> makeTransformator(real delta, real dx, real dy, real dz);
	virtual ~Transformator() {}

protected:
	Transformator(){}

public:
	virtual void transformWorldToGrid(Triangle &value) const = 0;
	virtual void transformWorldToGrid(TriangularMesh &geom) const = 0;
	virtual void transformWorldToGrid(Vertex &value) const = 0;

    virtual void transformGridToWorld(Triangle &t) const = 0;
	virtual void transformGridToWorld(Vertex &value) const = 0;
	
	virtual void transformGridToWorld(BoundingBox<real> &box) const = 0;
	virtual void transformWorldToGrid(BoundingBox<real> &box) const = 0;

};


#endif
