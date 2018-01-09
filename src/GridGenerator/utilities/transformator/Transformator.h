#ifndef Transformator_h
#define Transformator_h

#include <memory>
#include "GridGenerator/global.h"


template<typename T>
class BoundingBox;
struct Triangle;
struct Geometry;
struct Vertex;


class Transformator
{
public:
    static VF_PUBLIC std::shared_ptr<Transformator> makeTransformator(doubflo delta, doubflo dx, doubflo dy, doubflo dz);
	virtual ~Transformator() {}

protected:
	Transformator(){}

public:
	virtual void transformWorldToGrid(Triangle &value) const = 0;
	virtual void transformWorldToGrid(Geometry &geom) const = 0;
	virtual void transformWorldToGrid(Vertex &value) const = 0;

    virtual void transformGridToWorld(Triangle &t) const = 0;
	virtual void transformGridToWorld(Vertex &value) const = 0;
	
	virtual void transformGridToWorld(BoundingBox<doubflo> &box) const = 0;
	virtual void transformWorldToGrid(BoundingBox<doubflo> &box) const = 0;

};


#endif
