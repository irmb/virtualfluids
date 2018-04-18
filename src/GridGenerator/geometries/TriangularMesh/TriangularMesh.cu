#include "TriangularMesh.h"

#include <GridGenerator/io/STLReaderWriter/STLReader.h>
#include <GridGenerator/utilities/Transformator/TransformatorImp.h>
#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>
#include <time.h>

#include <utilities/logger/Logger.h>

#include "Timer/Timer.h"


TriangularMesh* TriangularMesh::make(const std::string& fileName, DiscretizationMethod discretizationMethod)
{
    TriangularMesh* triangularMesh = new TriangularMesh(fileName);
    for (uint i = 0; i < triangularMesh->triangleVec.size(); i++) {
        triangularMesh->minmax.setMinMax(triangularMesh->triangleVec[i]);
    }
    triangularMesh->discretizationMethod = discretizationMethod;
    return triangularMesh;
}

TriangularMesh::TriangularMesh(const std::string& input, const BoundingBox& box)
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
	this->minmax = BoundingBox(geo.minmax);
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

    Timer t = Timer::makeStart();

	TriangleNeighborFinder finder(triangles, size);
	finder.fillWithNeighborAngles(this);

    t.end();

    *logging::out << logging::Logger::INTERMEDIATE << "time finding neighbors: " << t.getTimeInSeconds() << "s\n";
}

void TriangularMesh::setTriangles(std::vector<Triangle> triangles)
{
	this->triangleVec = triangles;
	initalizeDataFromTriangles();
}

void TriangularMesh::setMinMax(BoundingBox minmax)
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


HOST DiscretizationMethod TriangularMesh::getDiscretizationMethod() const
{
    return this->discretizationMethod;
}