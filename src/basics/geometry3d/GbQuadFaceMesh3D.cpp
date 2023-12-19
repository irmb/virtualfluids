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
#include <geometry3d/GbQuadFaceMesh3D.h>

#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbHalfSpace3D.h>

using namespace std;

GbQuadFaceMesh3D::GbQuadFaceMesh3D() : GbObject3D()
{
    this->name       = "new GbMesh";
    this->nodes      = new vector<Vertex>;
    this->quads      = new vector<QuadFace>;
    this->consistent = false;
}

GbQuadFaceMesh3D::GbQuadFaceMesh3D(string name, vector<Vertex> *nodes, vector<QuadFace> *quads) : GbObject3D()
{
    if (name.size() == 0)
        throw UbException(UB_EXARGS, "no name specified");
    if (!nodes)
        throw UbException(UB_EXARGS, "no nodes specified");
    if (!quads)
        throw UbException(UB_EXARGS, "no quads specified");

    this->name       = name;
    this->nodes      = nodes;
    this->quads      = quads;
    this->consistent = false;
}
/*=============================================================================================*/

GbQuadFaceMesh3D::~GbQuadFaceMesh3D()
{
    if (nodes) {
        //    for(unsigned u=0; u<nodes->size(); u++) delete (*nodes)[u];
        delete nodes;
    }
    if (quads) {
        delete quads;
    }
}
/*======================================================================*/

void GbQuadFaceMesh3D::init()
{
    nodes      = NULL;
    quads      = NULL;
    x1min      = 0.0;
    x1max      = 0.0;
    x2min      = 0.0;
    x2max      = 0.0;
    x3min      = 0.0;
    x3max      = 0.0;
    consistent = false;
}
/**
 * Returns a string representation of this triangular mesh.
 * @return a string representation of this triangular mesh
 */
string GbQuadFaceMesh3D::toString()
{
    stringstream ss;
    ss << "GbQuadFaceMesh3D[";
    ss << (int)this->quads->size() << "-Quadangles, " << (int)this->nodes->size() << "-Nodes, " << endl;
    // ss<<"\""<<this->name<<", Area=sollt mal berechnet werden ;-)"<<"\"";
    // ss<<", x1min="<<this->x1min;
    // ss<<", x1max="<<this->x1max;
    // ss<<", x2min="<<this->x2min;
    // ss<<", x2max="<<this->x2max;
    // ss<<", x3min="<<this->x3min;
    // ss<<", x3max="<<this->x3max;
    ss << "]";
    return (ss.str());
}
/**
 * Returns the name of this triangular mesh.
 * @return the name of this triangular mesh
 */
string GbQuadFaceMesh3D::getName() { return (this->name); }

/**
 * Returns the nodes of this triangular mesh.
 * @return the nodes of this triangular mesh
 */
vector<GbQuadFaceMesh3D::Vertex> *GbQuadFaceMesh3D::getNodes() { return (this->nodes); }
/**
 * Returns the quads of this triangular mesh.
 * @return the quads of this triangular mesh
 */
vector<GbQuadFaceMesh3D::QuadFace> *GbQuadFaceMesh3D::getQuads() { return (this->quads); }
/**
 * Returns the center x1 coordinate of this triangular mesh.
 * @return the center x1 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX1Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (0.5 * (this->x1min + this->x1max));
}
/**
 * Returns the center x2 coordinate of this triangular mesh.
 * @return the center x2 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX2Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (0.5 * (this->x2min + this->x2max));
}
/**
 * Returns the center x3 coordinate of this triangular mesh.
 * @return the center x3 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX3Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (0.5 * (this->x3min + this->x3max));
}

/**
 * Returns the minimum x1 coordinate of this triangular mesh.
 * @return the minimum x1 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX1Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1min);
}
/**
 * Returns the maximum x1 coordinate of this triangular mesh.
 * @return the maximum x1 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX1Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1max);
}
/**
 * Returns the minimum x2 coordinate of this triangular mesh.
 * @return the minimum x2 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX2Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2min);
}
/**
 * Returns the maximum x2 coordinate of this triangular mesh.
 * @return the maximum x2 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX2Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2max);
}
/**
 * Returns the minimum x3 coordinate of this triangular mesh.
 * @return the minimum x3 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX3Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3min);
}
/**
 * Returns the maximum x3 coordinate of this triangular mesh.
 * @return the maximum x3 coordinate of this triangular mesh
 */
double GbQuadFaceMesh3D::getX3Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3max);
}

void GbQuadFaceMesh3D::calculateValues()
{
    double x1, x2, x3;

    this->x1min = (*this->nodes)[0].x;
    this->x1max = (*this->nodes)[0].x;
    this->x2min = (*this->nodes)[0].y;
    this->x2max = (*this->nodes)[0].y;
    this->x3min = (*this->nodes)[0].z;
    this->x3max = (*this->nodes)[0].z;

    for (int i = 1; i < (int)this->nodes->size(); i++) {
        x1 = (*this->nodes)[i].x;
        x2 = (*this->nodes)[i].y;
        x3 = (*this->nodes)[i].z;
        if (x1 < this->x1min)
            this->x1min = x1;
        if (x1 > this->x1max)
            this->x1max = x1;
        if (x2 < this->x2min)
            this->x2min = x2;
        if (x2 > this->x2max)
            this->x2max = x2;
        if (x3 < this->x3min)
            this->x3min = x3;
        if (x3 > this->x3max)
            this->x3max = x3;
    }
    this->consistent = true;
}

/*======================================================================*/
vector<GbTriangle3D *> GbQuadFaceMesh3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> triangles(0);
    return triangles;
    // throw UbException(__FILE__, __LINE__, "GbQuadFaceMesh3D::getSurfaceTriangelSet - not implemented");
}
// vector<GbQuad3D*> GbQuadFaceMesh3D::getSurfaceQuadSet()
//{
//   throw UbException(__FILE__, __LINE__, "GbQuadFaceMesh3D::getSurfaceQuadSet - not implemented");
//   //vector<GbQuadangle3D*> tris;
//   //GbQuadangle3D* quad;
//   //GbPoint3D* p1;
//   //GbPoint3D* p2;
//   //GbPoint3D* p3;
//   //int size = (int)this->quads->size();
//   //for(int u=0; u<size;u++)
//   //{
//   //   quad = (*this->quads)[u];
//   //   p1 = new GbPoint3D(quad->getPoint1());
//   //   p2 = new GbPoint3D(quad->getPoint2());
//   //   p3 = new GbPoint3D(quad->getPoint3());
//   //   tris.push_back(new GbQuadangle3D(p1, p2, p3));
//   //}
//   //return tris;
//}
/*======================================================================*/
/*
 * Function to determine if the point is inside the polyhedron defined as a 3D object
 * using the Halfspace algorithm
 * @param xp the x-coordinate of the point
 * @param yp the y-coordinate of the point
 * @param zp the z-coordinate of the point
 * @return true if point is inside else return false
 */
bool GbQuadFaceMesh3D::isPointInObject3DHalfSpace(const double & /*xp*/, const double & /*yp*/, const double & /*zp*/)
{
    throw UbException(UB_EXARGS, "not implemented");
    // vector<GbQuadangle3D*> *Quadangles = this->quads;
    // int Quadanglesize = (int)Quadangles->size();
    // GbPoint3D Point(xp,yp,zp);
    // for (int i=0; i<Quadanglesize; i++)
    //{
    //   GbPoint3D* point1 = (*Quadangles)[i]->getPoint1();
    //   GbPoint3D* point2 = (*Quadangles)[i]->getPoint2();
    //   GbPoint3D* point3 = (*Quadangles)[i]->getPoint3();

    //   GbHalfSpace3D halfspace(point1, point2, point3);
    //   if (halfspace.ptInside(&Point)) return false;
    //}
    // return true;
}
/*======================================================================*/
/*======================================================================*/
bool GbQuadFaceMesh3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3)
{

    double xmin = this->getX1Minimum();
    double xmax = this->getX1Maximum();
    double ymin = this->getX2Minimum();
    double ymax = this->getX2Maximum();
    double zmin = this->getX3Minimum();
    double zmax = this->getX3Maximum();
    double dX   = (xmax - xmin) / 100.;
    double dY   = (ymax - ymin) / 100.;
    double dZ   = (zmax - zmin) / 100.;
    GbCuboid3D boundingCube(xmin - dX, ymin - dY, zmin - dZ, xmax + dX, ymax + dY, zmax + dZ);
    if (!boundingCube.isPointInGbObject3D(x1, x2, x3)) {
        boundingCube.finalize();
        return false;
    }
    boundingCube.finalize();

    // Halfspace algorithm, Area of spherical polygons algorithm or Ray crossing algorithm
    GbVector3D bMin(boundingCube.getPoint1());
    GbVector3D bMax(boundingCube.getPoint2());
    bMin = bMax.Subtract(bMin);
    // int radius = (int)bMin.Length();

    // if(((GbQuadFaceMesh3D*)this->geoObject3D)->isPointInObject3DHalfSpace(x1,x2,x3) )
    // if(((GbQuadFaceMesh3D*)this->geoObject3D)->isPointInObject3Darea(x11,x12,x13,numQuadangles))
    // if(this->isPointInObject3DRayCrossing(x1,x2,x3,radius,(int)this->nodes->size(),(int)this->quads->size()))
    //   return true;
    // else
    return false;
}
/*======================================================================*/
bool GbQuadFaceMesh3D::isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/,
                                           bool & /*pointIsOnBoundary*/)
{
    throw UbException(UB_EXARGS, "not implemented");
}
/*======================================================================*/
GbLine3D *GbQuadFaceMesh3D::createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
}

//! \}
