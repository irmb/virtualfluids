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
    this->minmax = BoundingBox::makeInvalidMinMaxBox();

    this->triangleVec = STLReader::readSTL(inputPath);
    initalizeDataFromTriangles();
    this->findNeighbors();
}


TriangularMesh::TriangularMesh()
{

}

TriangularMesh::~TriangularMesh()
{

}

uint TriangularMesh::getNumberOfTriangles() const
{
    return triangleVec.size();
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
	this->size = long(triangleVec.size());

    for (uint i = 0; i < this->size; i++) {
        this->minmax.setMinMax(this->triangleVec[i]);
    }
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


bool intersectPlane(const Vertex &normal, const Vertex &pointOnPlane, const Vertex &originLine, const Vertex &directionLine, Vertex &intersectionPoint)
{
    // assuming vectors are all normalized
    real denom = normal * directionLine;
    if (denom > 1e-6) {
        Vertex p0l0 = pointOnPlane - originLine;
        real distance = p0l0 * normal / denom;
        intersectionPoint = originLine + directionLine * distance;
        return (distance >= 0);
    }

    return false;
}

void TriangularMesh::scale(double offset)
{
    std::vector<Triangle> triangles = this->triangleVec;

    TriangleNeighborFinder finder(this->triangles, this->size);

    auto trianglesPerVertex = finder.getTrianglesPerVertex();
    auto averrageNormals = getAverrageNormalsPerVertex(trianglesPerVertex);


    for (int vertexID = 0; vertexID < this->getNumberOfTriangles() * 3; vertexID++)
    {
        int coordinatedID = finder.sortedToTriangles[vertexID][IDS::coordinateID];
        Vertex averrageNormal = averrageNormals[coordinatedID];

        //auto triangleIds = finder.getTriangleIDsWithCommonVertex(vertexID);
        //for(auto index : triangleIds)
        //{
        //    auto triangle = origin->triangleVec[index];
        //    Vertex intersection;
        //    real d = 10;
        //    triangle.normal.normalize();
        //    Vertex p = triangle.v2 + triangle.normal * d;
        //    // p.normalize();
        //    Vertex lineOrigin = origin->triangleVec[triangleID].get(vertexTriangleID);
        //    // lineOrigin.normalize();
        //    //averrageNormal.normalize();
        //    //triangle.print();

        //    //printf("average: \n");
        //    //averrageNormal.print();

        //    //printf("line origin: \n");
        //    //lineOrigin.print();

        //    //printf("plane normal: \n");
        //    //triangle.normal.print();

        //    //printf("plane point: \n");
        //    //p.print();

        //    bool b = intersectPlane(triangle.normal, p, lineOrigin, averrageNormal, intersection);
        //    //intersection.print();
        //    //printf("\n");
        //}
        //printf("\n\n");

        averrageNormal.normalize();
        const int triangleID = vertexID / 3;
        const int vertexTriangleID = vertexID % 3;

        Vertex intersection;
        Vertex p = this->triangleVec[triangleID].v1 + this->triangleVec[triangleID].normal * offset;
        Vertex lineOrigin = this->triangleVec[triangleID].get(vertexTriangleID);
        bool b = intersectPlane(this->triangleVec[triangleID].normal, p, lineOrigin, averrageNormal, intersection);
        triangles[triangleID].set(vertexTriangleID, intersection);
        triangles[triangleID].calcNormal();

        triangles[triangleID].set(vertexTriangleID, lineOrigin + averrageNormal * offset);
    }

    //STLWriter::writeSTL(triangles, "D:/GRIDGENERATION/STL/FICKEN_ASCII.stl");
    this->setTriangles(triangles);
}

std::vector<Vertex> TriangularMesh::getAverrageNormalsPerVertex(std::vector<std::vector<Triangle>> trianglesPerVertex)
{
    std::vector<Vertex> averrageNormals;
    for (auto triangles : trianglesPerVertex)
    {
        eliminateTriangleswithIdenticialNormal(triangles);

        Vertex sumNormal;
        for (auto triangle : triangles)
        {
            triangle.calcNormal();
            sumNormal = sumNormal + triangle.normal;
        }
        real magnitude = sumNormal.getMagnitude();
        averrageNormals.push_back(Vertex(sumNormal / magnitude));
    }
    return averrageNormals;
}

void TriangularMesh::eliminateTriangleswithIdenticialNormal(std::vector<Triangle> &triangles)
{
    for (int i = 0; i < triangles.size() - 1; i++) {
        for (int j = i + 1; j < triangles.size(); j++) {
            if (triangles[i].normal == triangles[j].normal)
                triangles.erase(triangles.begin() + i);
        }
    }
}
