#include "TriangularMesh.h"

#include <GridGenerator/io/STLReaderWriter/STLReader.h>
#include <GridGenerator/utilities/Transformator/TransformatorImp.h>
#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>
#include <time.h>

#include <utilities/logger/Logger.h>

#include "Serialization/GeometryMemento.h"
#include <GridGenerator/geometries/Triangle/Serialization/TriangleMemento.h>
#include "Timer/Timer.h"


TriangularMesh* TriangularMesh::make(const std::string& fileName)
{
    TriangularMesh* triangularMesh = new TriangularMesh(fileName);
    for (uint i = 0; i < triangularMesh->triangleVec.size(); i++) {
        triangularMesh->minmax.setMinMax(triangularMesh->triangleVec[i]);
    }
    return triangularMesh;
}

TriangularMesh::TriangularMesh(const std::string& input, const BoundingBox<int>& box)
{
	this->triangleVec = STLReader::readSTL(box, input);
	initalizeDataFromTriangles();
	this->findNeighbors();
}

TriangularMesh::TriangularMesh(const std::string& inputPath)
{
    this->triangleVec = STLReader::readSTL(inputPath);
    initalizeDataFromTriangles();
    this->findNeighbors();
}

TriangularMesh::TriangularMesh(const TriangularMesh& geo)
{
	this->triangleVec = geo.triangleVec;
    this->triangles = geo.triangles;
    this->size = geo.size;
	this->minmax = BoundingBox<real>(geo.minmax);
}

TriangularMesh::TriangularMesh()
{

}

TriangularMesh::~TriangularMesh()
{

}


void TriangularMesh::findNeighbors()
{
	*logging::out << logging::Logger::INTERMEDIATE << "start finding neighbors ...\n";

	const clock_t begin = clock();
    Timer t = Timer::begin();

	TriangleNeighborFinder finder(triangles, size);
	finder.fillWithNeighborAngles(this);

	const clock_t end = clock();
    t.end();

	const real time = real(end - begin) / CLOCKS_PER_SEC;
    *logging::out << logging::Logger::INTERMEDIATE << "time finding neighbors: " << time << "s\n";
    *logging::out << logging::Logger::INTERMEDIATE << "time finding neighbors: " << t.getTimeInSeconds() << "s\n";
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


