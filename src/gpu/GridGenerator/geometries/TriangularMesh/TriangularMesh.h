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
//! \file TriangularMesh.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef TriangularMesh_h
#define TriangularMesh_h

#include <vector>
#include <string>
#include <memory>
#include "global.h"

#include "geometries/Triangle/Triangle.h"
#include "geometries/BoundingBox/BoundingBox.h"

#include "geometries/Object.h"

class GeometryMemento;
class GbTriFaceMesh3D;

enum class DiscretizationMethod { RAYCASTING, POINT_IN_OBJECT, POINT_UNDER_TRIANGLE };


class TriangularMesh : public Object
{
public:

    static TriangularMesh* make(const std::string& fileName, const std::vector<uint> ignorePatches = std::vector<uint>());
    TriangularMesh();
    TriangularMesh(const std::string& inputPath, const std::vector<uint> ignorePatches = std::vector<uint>());
    TriangularMesh(const std::string& inputPath, const BoundingBox &box);
    ~TriangularMesh() override = default;

    uint getNumberOfTriangles() const;

    void setTriangles(std::vector<Triangle> triangles);
    void setMinMax(BoundingBox minmax);

    std::vector<Triangle> triangleVec;
    Triangle *triangles = nullptr;
    long size = 0;
    BoundingBox minmax;

    SPtr<GbTriFaceMesh3D> VF_GbTriFaceMesh3D;

    bool operator==(const TriangularMesh &geometry) const;

    void findNeighbors();

    GbTriFaceMesh3D* getGbTriFaceMesh3D() const;

    void generateGbTriFaceMesh3D();

private:

    void initalizeDataFromTriangles();

    static std::vector<Vertex> getAverrageNormalsPerVertex(std::vector<std::vector<Triangle> > trianglesPerVertex);
    static void eliminateTriangleswithIdenticialNormal(std::vector<Triangle> &triangles);

public:
    SPtr<Object> clone() const override;
    double getX1Centroid() const override { throw "Not implemented in TriangularMesh"; }
    double getX1Minimum() const override { return minmax.minX; }
    double getX1Maximum() const override { return minmax.maxX; }
    double getX2Centroid() const override { throw "Not implemented in TriangularMesh"; }
    double getX2Minimum() const override { return minmax.minY; }
    double getX2Maximum() const override { return minmax.maxY; }
    double getX3Centroid() const override { throw "Not implemented in TriangularMesh"; }
    double getX3Minimum() const override { return minmax.minZ; }
    double getX3Maximum() const override { return minmax.maxZ; }
    void changeSizeByDelta(double delta) override;
    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset,
                         const double& maxOffset) override
    {
        (void)x1;
        (void)x2;
        (void)x3;
        (void)minOffset;
        (void)maxOffset;
        return false;
    }

    void findInnerNodes(SPtr<GridImp> grid) override;
};



#endif
