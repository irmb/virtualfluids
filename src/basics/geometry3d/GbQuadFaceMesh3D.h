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
#ifndef GBQUADFACEMESH3D_H
#define GBQUADFACEMESH3D_H

#include <iostream>
#include <sstream>

#include <basics/utilities/UbException.h>
#include <geometry3d/GbObject3D.h>

#include <PointerDefinitions.h>

class UbFileOutput;
class UbFileInput;
/*=========================================================================*/
/* GbQuadFaceMesh3D                                                                  */
/*                                                                         */
/**
 * This Class provides the triangular meshes.
 * Note, that up to now no methods for checking consistency are included.
 * in this context this class describes facettes from an 3D-object !!!
 */
class GbQuadFaceMesh3D : public GbObject3D
{
public:
    // nested class start
    class Vertex
    {
    public:
        Vertex() = default;
        Vertex(float x, float y, float z)
        {
            this->x = x;
            this->y = y;
            this->z = z;
        }
        float x, y, z;
    };

    class QuadFace
    {
    public:
        QuadFace() = default;
        QuadFace(int v1, int v2, int v3, int v4)
        {
            this->vertex1 = v1;
            this->vertex2 = v2;
            this->vertex3 = v3;
            this->vertex4 = v4;
        }

        int vertex1, vertex2, vertex3, vertex4;
    };
    // nested class end

public:
    GbQuadFaceMesh3D();
    GbQuadFaceMesh3D(std::string name, std::vector<Vertex> *nodes, std::vector<QuadFace> *quads);
    ~GbQuadFaceMesh3D() override;
    GbQuadFaceMesh3D *clone() override { throw UbException(UB_EXARGS, "clone() - not implemented"); }
    void finalize() override { throw UbException(UB_EXARGS, "finalize() - not implemented"); }

    std::string toString() override;
    std::string getName() override;
    std::vector<Vertex> *getNodes();
    std::vector<QuadFace> *getQuads();
    double getX1Centroid() override;
    double getX2Centroid() override;
    double getX3Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Minimum() override;
    double getX3Maximum() override;
    void calculateValues();

    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3) override;
    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary) override;

    bool isPointInObject3DHalfSpace(const double &xp, const double &yp,
                                    const double &zp); // based on Halfspace algorithm
    // bool isPointInObject3DSpherical(const double& xp, const double& yp, const double& zp, int numQuads);    //based
    // on Spherical polygon area method bool isPointInObject3DRayCrossing(const double& xp, const double& yp, const
    // double& zp, int radius, int numVertices, int numQuads);  //based on Ray tracing algorithm

    // char SegPlaneInt(GbQuad3D *quad, GbVector3D  &PointQ, GbVector3D &PointR, GbVector3D &Point, int *m);
    // char SegQuadCross(GbQuad3D *quad, GbVector3D  &PointQ, GbVector3D &PointR);
    // till here !!!

    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;
    // virtual std::vector<GbQuad3D*> getSurfaceQuadSet();
    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;

    virtual void write(UbFileOutput * /*out*/) { std::cout << "GbQuadFaceMesh3D::write - sorry not implemented\n"; }
    virtual void read(UbFileInput * /*in*/) { std::cout << "GbQuadFaceMesh3D::read  - sorry not implemented\n"; }

    void writeAVSMesh(UbFileOutput *out, bool normals = false);

    /*======================================================================*/
    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere
private:
    void init();
    /*======================================================================*/
    std::string name;
    std::vector<Vertex> *nodes;
    std::vector<QuadFace> *quads;
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
    bool consistent;
};
/*=========================================================================*/

#endif

//! \}
