#include "TriangularMesh.h"

#include <GridGenerator/io/STLReaderWriter/STLReader.h>
#include <GridGenerator/utilities/Transformator/TransformatorImp.h>
#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>
#include <time.h>

#include <utilities/logger/Logger.h>

#include "Serialization/GeometryMemento.h"
#include <GridGenerator/geometries/Triangle/Serialization/TriangleMemento.h>

TriangularMesh::TriangularMesh(const std::string& input, const BoundingBox<int>& box, const Transformator* trafo)
{
	this->transformator = new TransformatorImp(*(TransformatorImp*)(trafo));

	this->triangleVec = STLReader::readSTL(box, input, *trafo);
	initalizeDataFromTriangles();
	this->findNeighbors();
}

TriangularMesh::TriangularMesh(const std::string& inputPath)
{
    this->transformator = new TransformatorImp();

    this->triangleVec = STLReader::readSTL(inputPath);
    initalizeDataFromTriangles();
    this->findNeighbors();
}

TriangularMesh::TriangularMesh(const TriangularMesh& geo)
{
	this->transformator = new TransformatorImp();
	*this->transformator = *geo.transformator;
	this->triangleVec = geo.triangleVec;
    this->triangles = geo.triangles;
    this->size = geo.size;
	this->minmax = BoundingBox<real>(geo.minmax);
}

TriangularMesh::TriangularMesh()
{
	this->transformator = new TransformatorImp();
}

TriangularMesh::~TriangularMesh()
{
	delete this->transformator;
}

void TriangularMesh::transformChannelGeometry(const real resolution)
{
	delete this->transformator;
	this->transformator = new TransformatorImp(resolution, -minmax.minX, -minmax.minY, -minmax.minZ);

	transformator->transformWorldToGrid(*this);

	initalizeDataFromTriangles();
	findNeighbors();
}

void TriangularMesh::findNeighbors()
{
	*logging::out << logging::Logger::INTERMEDIATE << "start finding neighbors ...\n";

	clock_t begin = clock();
	TriangleNeighborFinder finder(triangles, size);
	finder.fillWithNeighborAngles(this);
	clock_t end = clock();

	real time = real(end - begin) / CLOCKS_PER_SEC;
	*logging::out << logging::Logger::INTERMEDIATE << "time finding neighbors: " << SSTR(time) << "s\n";
}

Transformator* TriangularMesh::getTransformator()
{
	return transformator;
}

void TriangularMesh::setTriangles(std::vector<Triangle> triangles)
{
	this->triangleVec = triangles;
	initalizeDataFromTriangles();
}

void TriangularMesh::setMinMax(BoundingBox<real> minmax)
{
	this->minmax = minmax;
}

void TriangularMesh::initalizeDataFromTriangles()
{
	this->triangles = triangleVec.data();
	this->size = (int)triangleVec.size();
}

HOST bool TriangularMesh::operator==(const TriangularMesh &geometry) const
{
    if (!(minmax == geometry.minmax))
        return false;

    if (size != geometry.size)
        return false;

    for (int i = 0; i < size ; i++)
        if (!(triangleVec[i] == geometry.triangleVec[i]))
            return false;
    return true;
}


HOST GeometryMemento TriangularMesh::getState() const
{
    GeometryMemento memento;
    for (int i = 0; i < size ; i++)
        memento.triangles.push_back(triangleVec[i].getState());
    
    return memento;
}

HOST void TriangularMesh::setState(const GeometryMemento &memento)
{
    this->size = memento.triangles.size();
    this->triangleVec.resize(size);

    Triangle t;
    for (int i = 0; i < size; i++) {
        t.setState(memento.triangles[i]);
        triangleVec[i] = t;
    }
    this->initalizeDataFromTriangles();
}


