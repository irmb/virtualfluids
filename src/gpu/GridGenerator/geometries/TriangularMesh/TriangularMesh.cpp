//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TriangularMesh.cpp
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "TriangularMesh.h"

#include "Timer/Timer.h"

#include "basics/geometry3d/GbTriFaceMesh3D.h"

#include "geometries/TriangularMesh/triangleNeighborFinder/TriangleNeighborFinder.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "io/STLReaderWriter/STLWriter.h"
#include "io/STLReaderWriter/STLReader.h"

#include "grid/GridImp.h"
#include "grid/NodeValues.h"

#include "utilities/transformator/TransformatorImp.h"

using namespace vf::gpu;


TriangularMesh* TriangularMesh::make(const std::string& fileName, const std::vector<uint> ignorePatches)
{
    TriangularMesh* triangularMesh = new TriangularMesh(fileName, ignorePatches);
    return triangularMesh;
}

TriangularMesh::TriangularMesh(const std::string& input, const BoundingBox& box)
{
    this->triangleVec = STLReader::readSTL(box, input);
    initalizeDataFromTriangles();
    this->findNeighbors();
}

TriangularMesh::TriangularMesh(const std::string& inputPath, const std::vector<uint> ignorePatches)
{
    this->minmax = BoundingBox::makeInvalidMinMaxBox();

    this->triangleVec = STLReader::readSTL(inputPath, STLReader::ascii, ignorePatches);
    //this->triangleVec = STLReader::readSTL(inputPath);
    initalizeDataFromTriangles();
    this->findNeighbors();
}


TriangularMesh::TriangularMesh()
{
    this->minmax = BoundingBox::makeInvalidMinMaxBox();  // blame Lenz
}

SPtr<Object> TriangularMesh::clone() const
{
    auto mesh = std::make_shared<TriangularMesh>();
    mesh->setTriangles(this->triangleVec);
    return mesh;
}


uint TriangularMesh::getNumberOfTriangles() const
{
    return (uint)triangleVec.size();
}


void TriangularMesh::findNeighbors()
{
    VF_LOG_INFO("start finding neighbors ...");

    vf::basics::Timer t;
    t.start();

    TriangleNeighborFinder finder(triangles, size);
    finder.fillWithNeighborAngles(this);

    t.end();
    VF_LOG_INFO("time finding neighbors = {}", t.getTimeInSeconds());
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

    for (std::size_t i = 0; i < (size_t)this->size; i++) {
        this->minmax.setMinMax(this->triangleVec[i]);
    }
}

bool TriangularMesh::operator==(const TriangularMesh &geometry) const
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


GbTriFaceMesh3D* TriangularMesh::getGbTriFaceMesh3D() const
{
    return this->VF_GbTriFaceMesh3D.get();
}

void TriangularMesh::generateGbTriFaceMesh3D()
{
    if( this->VF_GbTriFaceMesh3D ) return;

    VF_LOG_INFO("Start generating GbTriFaceMesh3D");

    std::vector<GbTriFaceMesh3D::Vertex>  *gbVertices = new std::vector<GbTriFaceMesh3D::Vertex>(this->triangleVec.size() * 3);
    std::vector<GbTriFaceMesh3D::TriFace> *gbTriangles = new std::vector<GbTriFaceMesh3D::TriFace>(this->triangleVec.size());

    for (int i = 0; i < (int)this->triangleVec.size(); i++)
    {
        (*gbVertices)[i * 3] = GbTriFaceMesh3D::Vertex(triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
        (*gbVertices)[i * 3 + 1] = GbTriFaceMesh3D::Vertex(triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
        (*gbVertices)[i * 3 + 2] = GbTriFaceMesh3D::Vertex(triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);

        (*gbTriangles)[i] = GbTriFaceMesh3D::TriFace(i * 3, i * 3 + 1, i * 3 + 2);
    }

    this->VF_GbTriFaceMesh3D = std::make_shared<GbTriFaceMesh3D>( "stl", gbVertices, gbTriangles, GbTriFaceMesh3D::KDTREE_SAHPLIT, false );

    VF_LOG_INFO("Done generating GbTriFaceMesh3D");
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

void TriangularMesh::changeSizeByDelta(double offset)
{
    std::vector<Triangle> triangles = this->triangleVec;

    TriangleNeighborFinder finder(this->triangles, this->size);

    auto trianglesPerVertex = finder.getTrianglesPerVertex();
    auto averrageNormals = getAverrageNormalsPerVertex(trianglesPerVertex);


    for (uint vertexID = 0; vertexID < this->getNumberOfTriangles() * 3; vertexID++)
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
        const int triangleID = (int)vertexID / 3;
        const int vertexTriangleID = (int)vertexID % 3;

        Vertex intersection;
        // Vertex p = this->triangleVec[triangleID].v1 + this->triangleVec[triangleID].normal * offset; // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/85
        Vertex lineOrigin = this->triangleVec[triangleID].get(vertexTriangleID);
        // bool b = intersectPlane(this->triangleVec[triangleID].normal, p, lineOrigin, averrageNormal, intersection);
        triangles[triangleID].set(vertexTriangleID, intersection);
        triangles[triangleID].calcNormal();

        triangles[triangleID].set(vertexTriangleID, lineOrigin + averrageNormal * offset);
    }

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
    for (std::size_t i = 0; i < triangles.size() - 1; i++) {
        for (std::size_t j = i + 1; j < triangles.size(); j++) {
            if (triangles[i].normal == triangles[j].normal)
                triangles.erase(triangles.begin() + i);
        }
    }
}

void TriangularMesh::findInnerNodes(SPtr<GridImp> grid)
{
    grid->getTriangularMeshDiscretizationStrategy()->discretize(this, grid.get(), FLUID, INVALID_OUT_OF_GRID);
}
