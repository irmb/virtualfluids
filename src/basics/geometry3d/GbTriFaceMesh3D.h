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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Konstantin Kutscher, Soeren Textor, Sebastian Geller
//=======================================================================================
#ifndef GBTRIFACEMESH3D_H
#define GBTRIFACEMESH3D_H

#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/Vector3D.h>

#include <geometry3d/GbPoint3D.h>

#include <PointerDefinitions.h>

#include "basics/constants/NumericConstants.h"

namespace Kd
{
template <typename T>
class Tree;
template <typename T>
class SplitAlgorithm;
template <typename T>
class RayIntersectionHandler;
} // namespace Kd

class WbWriter;

/*=========================================================================*/
/* GbTriFaceMesh3D                                                                  */
/*                                                                         */
/**
 * This Class provides the triangular meshes.
 * Note, that up to now no methods for checking consistency are included.
 * in this context this class describes facettes from an 3D-object !!!
 */
class GbTriFaceMesh3D : public GbObject3D
{
public:
    // nested class start
    class Vertex
    {
    public:
        Vertex() = default;
        Vertex(const float &x, const float &y, const float &z) : x(x), y(y), z(z) {}
        Vertex(Vertex *vert)
        {
            this->x = vert->x;
            this->y = vert->y;
            this->z = vert->z;
        }
        float operator[](const int &i) const
        {
            if (i == 0)
                return x;
            else if (i == 1)
                return y;
            else if (i == 2)
                return z;

            throw UbException(UB_EXARGS, "i not in [0;2]");
        }
        float &operator[](const int &i)
        {
            if (i == 0)
                return x;
            else if (i == 1)
                return y;
            else if (i == 2)
                return z;

            throw UbException(UB_EXARGS, "not in [0;2]");
        }
        bool operator==(const Vertex &rhs)
        {
            return (fabs(x - rhs.x) < 1.E-8 && fabs(y - rhs.y) < 1.E-8 && fabs(z - rhs.z) < 1.E-8);
        }
        friend inline bool operator<(const Vertex &lhsVert, const Vertex &rhsVert)
        {
            if (lhsVert.x < rhsVert.x)
                return true;
            if (lhsVert.x > rhsVert.x)
                return false;
            if (lhsVert.y < rhsVert.y)
                return true;
            if (lhsVert.y > rhsVert.y)
                return false;
            if (lhsVert.z < rhsVert.z)
                return true;

            return false;
        }
        friend std::ostream &operator<<(std::ostream &os, const Vertex &node)
        {
            return os << node.x << "," << node.y << "," << node.z;
        }
        Vertex *clone() { return (new Vertex(this)); }

    public:
        float x{ 0.0 }, y{ 0.0 }, z{ 0.0 };
    };
    //////////////////////////////////////////////////////////////////////////
    class TriFace
    {
    public:
        TriFace()

            = default;
        TriFace(const int &v1, const int &v2, const int &v3) : v1(v1), v2(v2), v3(v3) {}

        const int &getIndexVertex1() const { return v1; }
        const int &getIndexVertex2() const { return v2; }
        const int &getIndexVertex3() const { return v3; }

        Vertex &getNode(const int &i, std::vector<Vertex> &nodes)
        {
            if (i == 0)
                return nodes[v1];
            if (i == 1)
                return nodes[v2];
            if (i == 2)
                return nodes[v3];
            throw UbException(UB_EXARGS, "invalid i - not in range [0;2]");
        }
        void setNode(const int &i, const int &index)
        {
            if (i == 0)
                v1 = index;
            else if (i == 1)
                v2 = index;
            else if (i == 2)
                v3 = index;
            else
                throw UbException(UB_EXARGS, "invalid i - not in range [0;2]");
        }

        int operator[](int index)
        {
            if (index == 0)
                return v1;
            if (index == 1)
                return v2;
            if (index == 2)
                return v3;
            throw UbException(UB_EXARGS, "invalid i - not in range [0;2]");
        }

        float &getV1x(std::vector<Vertex> &nodes) { return nodes[v1].x; }
        float &getV1y(std::vector<Vertex> &nodes) { return nodes[v1].y; }
        float &getV1z(std::vector<Vertex> &nodes) { return nodes[v1].z; }

        float &getV2x(std::vector<Vertex> &nodes) { return nodes[v2].x; }
        float &getV2y(std::vector<Vertex> &nodes) { return nodes[v2].y; }
        float &getV2z(std::vector<Vertex> &nodes) { return nodes[v2].z; }

        float &getV3x(std::vector<Vertex> &nodes) { return nodes[v3].x; }
        float &getV3y(std::vector<Vertex> &nodes) { return nodes[v3].y; }
        float &getV3z(std::vector<Vertex> &nodes) { return nodes[v3].z; }

        float getMinX(std::vector<Vertex> &nodes) { return (float)UbMath::min(nodes[v1].x, nodes[v2].x, nodes[v3].x); }
        float getMinY(std::vector<Vertex> &nodes) { return (float)UbMath::min(nodes[v1].y, nodes[v2].y, nodes[v3].y); }
        float getMinZ(std::vector<Vertex> &nodes) { return (float)UbMath::min(nodes[v1].z, nodes[v2].z, nodes[v3].z); }

        float getMaxX(std::vector<Vertex> &nodes) { return (float)UbMath::max(nodes[v1].x, nodes[v2].x, nodes[v3].x); }
        float getMaxY(std::vector<Vertex> &nodes) { return (float)UbMath::max(nodes[v1].y, nodes[v2].y, nodes[v3].y); }
        float getMaxZ(std::vector<Vertex> &nodes) { return (float)UbMath::max(nodes[v1].z, nodes[v2].z, nodes[v3].z); }

        float getX1Centroid(std::vector<Vertex> &nodes)
        {
            return (float)vf::basics::constant::c1o3 * (getV1x(nodes) + getV2x(nodes) + getV3x(nodes));
        }
        float getX2Centroid(std::vector<Vertex> &nodes)
        {
            return (float)vf::basics::constant::c1o3 * (getV1y(nodes) + getV2y(nodes) + getV3y(nodes));
        }
        float getX3Centroid(std::vector<Vertex> &nodes)
        {
            return (float)vf::basics::constant::c1o3 * (getV1z(nodes) + getV2z(nodes) + getV3z(nodes));
        }

        double calculateDistanceToPoint3D(const double &x1, const double &x2, const double &x3,
                                          std::vector<Vertex> &nodes);

        double getArea(std::vector<Vertex> &nodes)
        {
            // GbVector3D A(nodes[v1].x, nodes[v1].y, nodes[v1].z);
            // GbVector3D B(nodes[v2].x, nodes[v2].y, nodes[v2].z);
            // GbVector3D C(nodes[v3].x, nodes[v3].y, nodes[v3].z);
            // GbVector3D AB = B-A;
            // GbVector3D AC = C-A;
            // GbVector3D N = AB.Cross(AC);
            // return 0.5*N.Length();
            Vector3D A(nodes[v1].x, nodes[v1].y, nodes[v1].z);
            Vector3D B(nodes[v2].x, nodes[v2].y, nodes[v2].z);
            Vector3D C(nodes[v3].x, nodes[v3].y, nodes[v3].z);
            Vector3D AB = B - A;
            Vector3D AC = C - A;
            Vector3D N  = AB.Cross(AC);
            return 0.5 * N.Length();
        }
        void calculateNormal(std::vector<Vertex> &nodes)
        {
            const float &v1x = nodes[v1].x;
            const float &v1y = nodes[v1].y;
            const float &v1z = nodes[v1].z;
            const float &v2x = nodes[v2].x;
            const float &v2y = nodes[v2].y;
            const float &v2z = nodes[v2].z;
            const float &v3x = nodes[v3].x;
            const float &v3y = nodes[v3].y;
            const float &v3z = nodes[v3].z;

            nx = (v3z - v1z) * (v2y - v1y) - (v2z - v1z) * (v3y - v1y);
            ny = (v2z - v1z) * (v3x - v1x) - (v2x - v1x) * (v3z - v1z);
            nz = (v2x - v1x) * (v3y - v1y) - (v2y - v1y) * (v3x - v1x);

            float length = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (length > 1.E-10) {
                length = 1.0f / length;
                nx *= length;
                ny *= length;
                nz *= length;
            } else {
                std::cerr << "GbTriFaceMesh3D::TriFace - calculateNormal: nx=ny=nz=0 -> kann nich sein "
                          << "(dreieck hat evtl knoten doppelt oder ist ne Linie)"
                          << "->removeRedunantNodes" << std::endl;
            }
        }

    public:
        int v1{ -1 }, v2{ -1 }, v3{ -1 };
        float nx{ 0.0 }, ny{ 0.0 }, nz{ 0.0 };
    };

public:
    enum KDTREE_SPLITAGORITHM { KDTREE_SAHPLIT, KDTREE_SPATIALSPLIT };

public:
    GbTriFaceMesh3D();
    GbTriFaceMesh3D(std::string name, std::vector<Vertex> *nodes, std::vector<TriFace> *triangles,
                    KDTREE_SPLITAGORITHM splitAlg = KDTREE_SAHPLIT, bool removeRedundantNodes = true);
    ~GbTriFaceMesh3D() override;

    GbTriFaceMesh3D *clone() override; // { throw UbException(UB_EXARGS,"not implemented"); }
    void finalize() override {}

    // void setRegardPointInPolyhedronTest(bool value) { this->regardPiO=value; }

    std::string toString() override;

    // std::string getName();
    std::vector<Vertex> *getNodes();
    std::vector<TriFace> *getTriangles();

    void setTransferViaFilename(bool transferViaFilename, std::string filename, double transX1, double transX2,
                                double transX3)
    {
        this->filename            = filename;
        this->transferViaFilename = transferViaFilename;
        this->transX1             = transX1;
        this->transX2             = transX2;
        this->transX3             = transX3;
    }
    void readMeshFromSTLFileASCII(std::string filename, bool removeRedundantNodes);
    void readMeshFromSTLFileBinary(std::string filename, bool removeRedundantNodes);

    double getX1Minimum() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x1min;
    }
    double getX1Maximum() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x1max;
    }
    double getX1Centroid() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x1center;
    }

    double getX2Minimum() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x2min;
    }
    double getX2Maximum() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x2max;
    }
    double getX2Centroid() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x2center;
    }

    double getX3Minimum() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x3min;
    }
    double getX3Centroid() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x3center;
    }
    double getX3Maximum() override
    {
        if (!this->consistent)
            this->calculateValues();
        return this->x3max;
    }

    void calculateValues();

    double getVolume();
    void deleteRedundantNodes();

    UbTupleDouble6 calculateMomentOfInertia(double rhoP);
    UbTupleDouble3 calculateCenterOfGravity();

    void setCenterCoordinates(const double &x1, const double &x2, const double &x3) override;

    void setCenterCoordinates(const UbTupleDouble3 & /*position*/) override;

    void scale(const double &sx1, const double &sx2, const double &sx3) override;
    void rotate(const double &alpha, const double &beta, const double &gamma) override;
    void rotateAroundPoint(const double &px1, const double &px2, const double &px3, const double &alpha,
                           const double &beta, const double &gamma);
    void translate(const double &x1, const double &x2, const double &x3) override;
    //void reflectAcrossXYLine(const double &alpha);

    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3) override;
    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3, int counter);
    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary) override;

    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles) override;

    std::vector<GbTriFaceMesh3D::TriFace *> getTrianglesForVertex(Vertex *vertex);

    void setKdTreeSplitAlgorithm(KDTREE_SPLITAGORITHM mode);
    KDTREE_SPLITAGORITHM getKdTreeSplitAlgorithm() { return this->kdtreeSplitAlg; }
    Kd::Tree<double> *getKdTree() { return this->kdTree; }

    virtual UbTuple<std::string, std::string> writeMesh(std::string filename, WbWriter *writer,
                                                        bool writeNormals                          = false,
                                                        std::vector<std::string> *datanames        = NULL,
                                                        std::vector<std::vector<double>> *nodedata = NULL);
    void writeMeshPly(const std::string &filename);

    /*======================================================================*/
    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere

    bool intersectLine(const double &p1_x1, const double &p1_x2, const double &p1_x3, const double &p2_x1,
                       const double &p2_x2, const double &p2_x3);

protected:
    KDTREE_SPLITAGORITHM kdtreeSplitAlg;
    void init();

    std::vector<Vertex> * nodes;
    std::vector<TriFace> * triangles;
    // for transfer
    std::string filename;
    bool transferViaFilename{ false };
    double transX1{ 0.0 };
    double transX2{ 0.0 };
    double transX3{ 0.0 };

    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
    double x1center;
    double x2center;
    double x3center;

    bool consistent{ false };

    bool buildVertTriRelationMap{ false };
    std::multimap<Vertex *, TriFace *> relationVertTris;

    Kd::Tree<double> *kdTree = nullptr;
};

#endif // GBTRIFACEMESH3D_H

//! \}
