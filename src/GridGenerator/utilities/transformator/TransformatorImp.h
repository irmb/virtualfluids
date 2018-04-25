#ifndef TransformatorImp_h
#define TransformatorImp_h


#include "GridGenerator/global.h"
#include <exception>
#include <sstream>
#include "Transformator.h"
#include "ArrowTransformator.h"

class BoundingBox;
struct Triangle;
class TriangularMesh;
struct Vertex;


class invalidDelta : public std::exception
{
	const char* what() const throw() {
		std::ostringstream getNr;
		getNr << "Delta cant be < Null. To enable no changes change delta to 1.0.";
		return getNr.str().c_str();
	}
};

class TransformatorImp
	: public Transformator, public ArrowTransformator
{
public:
	VF_PUBLIC TransformatorImp();
	VF_PUBLIC TransformatorImp(const TransformatorImp& trafo);
	VF_PUBLIC TransformatorImp(real delta, Vertex& translater);
	VF_PUBLIC TransformatorImp(real delta, real dx, real dy, real dz);
	VF_PUBLIC virtual ~TransformatorImp();
	
	VF_PUBLIC void transformWorldToGrid(Triangle &value) const;
	VF_PUBLIC void transformWorldToGrid(TriangularMesh &geom) const;
	VF_PUBLIC void transformWorldToGrid(Vertex &value) const;

    VF_PUBLIC void transformGridToWorld(Triangle &t) const;
	VF_PUBLIC void transformGridToWorld(Vertex &value) const;

	VF_PUBLIC void transformGridToWorld(BoundingBox &box) const;
	VF_PUBLIC void transformWorldToGrid(BoundingBox &box) const;

	VF_PUBLIC bool operator==(const TransformatorImp& trafo) const;

	VF_PUBLIC virtual void transformGridToWorld(std::shared_ptr<Arrow> arrow) const override;

private:
	real delta;
	std::shared_ptr<Vertex> translater;

	void scaleWorldToView(Vertex & v) const;
	void translateWorldToView(Vertex & v) const;

	void translateGridToWorld(Vertex & value) const;
	void scaleGridToWorld(Vertex & value) const;

	void verifyDelta(real delta) const;
};


#endif
