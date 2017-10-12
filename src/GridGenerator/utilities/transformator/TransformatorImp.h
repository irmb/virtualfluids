#ifndef TransformatorImp_h
#define TransformatorImp_h

#include "GridGenerator_EXPORT.h"
#include "GridGenerator/global.h"
#include <exception>

#include "Transformator.h"
#include "ArrowTransformator.h"

template<typename T>
class BoundingBox;
struct Triangle;
struct Geometry;
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
	GridGenerator_EXPORT TransformatorImp();
	GridGenerator_EXPORT TransformatorImp(const TransformatorImp& trafo);
	GridGenerator_EXPORT TransformatorImp(doubflo delta, Vertex& translater);
	GridGenerator_EXPORT TransformatorImp(doubflo delta, doubflo dx, doubflo dy, doubflo dz);
	GridGenerator_EXPORT virtual ~TransformatorImp();
	
	GridGenerator_EXPORT void transformWorldToGrid(Triangle &value) const;
	GridGenerator_EXPORT void transformWorldToGrid(Geometry &geom) const;
	GridGenerator_EXPORT void transformWorldToGrid(Vertex &value) const;

    GridGenerator_EXPORT void transformGridToWorld(Triangle &t) const;
	GridGenerator_EXPORT void transformGridToWorld(Vertex &value) const;

	GridGenerator_EXPORT void transformGridToWorld(BoundingBox<doubflo> &box) const;
	GridGenerator_EXPORT void transformWorldToGrid(BoundingBox<doubflo> &box) const;

	GridGenerator_EXPORT bool operator==(const TransformatorImp& trafo) const;

	GridGenerator_EXPORT virtual void transformGridToWorld(std::shared_ptr<Arrow> arrow) const override;

private:
	doubflo delta;
	std::shared_ptr<Vertex> translater;

	void scaleWorldToView(Vertex & v) const;
	void translateWorldToView(Vertex & v) const;

	void translateGridToWorld(Vertex & value) const;
	void scaleGridToWorld(Vertex & value) const;

	void verifyDelta(doubflo delta) const;
};


#endif
