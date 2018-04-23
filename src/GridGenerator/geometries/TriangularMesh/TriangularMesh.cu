#include "TriangularMesh.h"

#include <GridGenerator/io/STLReaderWriter/STLReader.h>
#include <GridGenerator/utilities/Transformator/TransformatorImp.h>
#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>
#include <time.h>

#include <utilities/logger/Logger.h>

#include "Timer/Timer.h"

#include "numerics/geometry3d/GbTriFaceMesh3D.h"


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
    this->minmax = BoundingBox::makeInvalidMinMaxBox();
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

HOSTDEVICE GbTriFaceMesh3D* TriangularMesh::getGbTriFaceMesh3D() const
{
    std::vector<GbTriFaceMesh3D::Vertex> *gbVertices = new std::vector<GbTriFaceMesh3D::Vertex>(this->triangleVec.size() * 3);
    std::vector<GbTriFaceMesh3D::TriFace> *gbTriangles = new std::vector<GbTriFaceMesh3D::TriFace>(this->triangleVec.size());
    for (int i = 0; i < this->triangleVec.size(); i++)
    {
        (*gbVertices)[i * 3] = GbTriFaceMesh3D::Vertex(triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
        (*gbVertices)[i * 3 + 1] = GbTriFaceMesh3D::Vertex(triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
        (*gbVertices)[i * 3 + 2] = GbTriFaceMesh3D::Vertex(triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);

        (*gbTriangles)[i] = GbTriFaceMesh3D::TriFace(i * 3, i * 3 + 1, i * 3 + 2);
    }

    GbTriFaceMesh3D* mesh = new GbTriFaceMesh3D("stl", gbVertices, gbTriangles);
    return mesh;
}