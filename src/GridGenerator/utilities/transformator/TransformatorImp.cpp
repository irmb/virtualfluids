#include "TransformatorImp.h"
#include <memory>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>
#include <GridGenerator/geometries/Triangle/Triangle.h>
#include <GridGenerator/geometries/TriangularMesh/TriangularMesh.h>
#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include <GridGenerator/geometries/Arrow/Arrow.h>

TransformatorImp::TransformatorImp() 
{
	this->translater = std::shared_ptr<Vertex>(new Vertex());
	this->delta = 1.0;
	this->translater->x = 0;
	this->translater->y = 0;
	this->translater->z = 0;
}

TransformatorImp::TransformatorImp(real delta, Vertex& translater) : delta(delta), translater(std::make_shared<Vertex>(translater))
{
	this->verifyDelta(delta);
}

TransformatorImp::TransformatorImp(real delta, real dx, real dy, real dz) : TransformatorImp(delta, Vertex(dx,dy,dz))
{

}

TransformatorImp::TransformatorImp(const TransformatorImp& trafo)
{
	this->delta = trafo.delta;
	this->translater = std::shared_ptr<Vertex>(new Vertex(*trafo.translater.get()));
}

TransformatorImp::~TransformatorImp()
{

}

void TransformatorImp::transformWorldToGrid(TriangularMesh &geom) const
{
	for (int i = 0; i < geom.size; i++)
		transformWorldToGrid(geom.triangleVec[i]);
}

void TransformatorImp::transformWorldToGrid(Triangle &value) const
{
	transformWorldToGrid(value.v1);
	transformWorldToGrid(value.v2);
	transformWorldToGrid(value.v3);
}

void TransformatorImp::transformGridToWorld(std::shared_ptr<Arrow> arrow) const
{
	transformGridToWorld(*arrow->getStart());
	transformGridToWorld(*arrow->getEnd());
}

void TransformatorImp::transformWorldToGrid(Vertex &v) const
{
	translateWorldToView(v);
	scaleWorldToView(v);
}


void TransformatorImp::translateWorldToView(Vertex& v) const
{
	v = *translater.get() + v;
}

void TransformatorImp::scaleWorldToView(Vertex& v) const
{
	v = v * (1.0f / this->delta);
}


void TransformatorImp::transformGridToWorld(Triangle & t) const
{
    transformGridToWorld(t.v1);
    transformGridToWorld(t.v2);
    transformGridToWorld(t.v3);
}

void TransformatorImp::transformGridToWorld(Vertex &value) const
{
	scaleGridToWorld(value);
	translateGridToWorld(value);
}

void TransformatorImp::scaleGridToWorld(Vertex & value) const
{
	value = value * this->delta;
}


void TransformatorImp::translateGridToWorld(Vertex & value) const
{
	value = value - *translater.get();
}


void TransformatorImp::transformGridToWorld(BoundingBox<real> &box) const
{
	//scale
	box.minX = (box.minX * this->delta);
	box.minY = (box.minY * this->delta);
	box.minZ = (box.minZ * this->delta);

	box.maxX = (box.maxX * this->delta);
	box.maxY = (box.maxY * this->delta);
	box.maxZ = (box.maxZ * this->delta);

	//translate
	box.minX = (box.minX - this->translater->x);
	box.minY = (box.minY - this->translater->y);
	box.minZ = (box.minZ - this->translater->z);

	box.maxX = (box.maxX - this->translater->x);
	box.maxY = (box.maxY - this->translater->y);
	box.maxZ = (box.maxZ - this->translater->z);
}

void TransformatorImp::transformWorldToGrid(BoundingBox<real> &box) const
{
	//translate
	box.minX += this->translater->x;
	box.minY += this->translater->y;
	box.minZ += this->translater->z;

	box.maxX += this->translater->x;
	box.maxY += this->translater->y;
	box.maxZ += this->translater->z;

	//scale
	box.minX *= (1.0f / this->delta);
	box.minY *= (1.0f / this->delta);
	box.minZ *= (1.0f / this->delta);

	box.maxX *= (1.0f / this->delta);
	box.maxY *= (1.0f / this->delta);
	box.maxZ *= (1.0f / this->delta);
}


void TransformatorImp::verifyDelta(real delta) const
{
	if (delta <= 0.0)
		throw invalidDelta();
}

bool TransformatorImp::operator==(const TransformatorImp& trafo) const
{
	return (this->delta == trafo.delta
		&& this->translater->x == trafo.translater->x
		&& this->translater->y == trafo.translater->y
		&& this->translater->z == trafo.translater->z);
}
