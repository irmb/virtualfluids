#ifndef TransformatorMocks_h
#define TransformatorMocks_h

#include "gmock/gmock.h"


#include "GridGenerator/global.h"

#include "Transformator.h"

#include <geometries/Vertex/Vertex.cuh>
#include <geometries/Triangle/Triangle.h>
#include <geometries/BoundingBox/BoundingBox.h>
#include <geometries/TriangularMesh/TriangularMesh.h>

class TransformatorStub : public Transformator
{
public:
	TransformatorStub() {}
	virtual ~TransformatorStub() {};

	virtual void transformWorldToView(Triangle &value) const override{}
	virtual void transformWorldToView(TriangularMesh &geom) const override {}
	virtual void transformWorldToView(Vertex &value) const override {}
	
	virtual void transformViewToWorld(BoundingBox<real> &box) const override {}
	virtual void transformWorldToView(BoundingBox<real> &box) const override {}

	virtual void transformViewToWorld(Vertex &value) const override { v = value; counter++; logString.append("transformViewToWorld "); }

	Vertex getParameterFromTransViewToWorld() { return v; }
	std::string getLogString() { return logString; }

private:
	mutable Vertex v;
	mutable int counter = 0;
	mutable std::string logString;
};

class TransformatorSpy : public TransformatorStub
{
public:
	virtual ~TransformatorSpy() {}
	MOCK_CONST_METHOD1(transformViewToWorld, void(Vertex& value));

};



#endif
