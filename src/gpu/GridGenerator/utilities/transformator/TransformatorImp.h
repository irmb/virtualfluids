#ifndef TransformatorImp_h
#define TransformatorImp_h

#include <exception>
#include <sstream>

#include "global.h"
#include "GridGenerator_export.h"

#include "utilities/transformator/Transformator.h"
#include "utilities/transformator/ArrowTransformator.h"

class BoundingBox;
struct Triangle;
class TriangularMesh;
struct Vertex;

class invalidDelta : public std::exception
{
    public:
    invalidDelta() : error_message ("Delta cant be < Null. To enable no changes change delta to 1.0.")
    {}

	const char* what() const noexcept
    {
	    return error_message.c_str();
	}

    private:
    std::string error_message;
};

class TransformatorImp
	: public Transformator, public ArrowTransformator
{
public:
	GRIDGENERATOR_EXPORT TransformatorImp();
	GRIDGENERATOR_EXPORT TransformatorImp(const TransformatorImp& trafo);
	GRIDGENERATOR_EXPORT TransformatorImp(real delta, const Vertex& translater);
	GRIDGENERATOR_EXPORT TransformatorImp(real delta, real dx, real dy, real dz);
	GRIDGENERATOR_EXPORT virtual ~TransformatorImp();
	
	GRIDGENERATOR_EXPORT void transformWorldToGrid(Triangle &value) const override;
	GRIDGENERATOR_EXPORT void transformWorldToGrid(TriangularMesh &geom) const override;
	GRIDGENERATOR_EXPORT void transformWorldToGrid(Vertex &value) const override;

    GRIDGENERATOR_EXPORT void transformGridToWorld(Triangle &t) const override;
	GRIDGENERATOR_EXPORT void transformGridToWorld(Vertex &value) const override;

	GRIDGENERATOR_EXPORT void transformGridToWorld(BoundingBox &box) const override;
	GRIDGENERATOR_EXPORT void transformWorldToGrid(BoundingBox &box) const override;

	GRIDGENERATOR_EXPORT bool operator==(const TransformatorImp& trafo) const;

	GRIDGENERATOR_EXPORT virtual void transformGridToWorld(std::shared_ptr<Arrow> arrow) const override;

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
