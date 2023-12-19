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
#include <geometry3d/GbTriangularMesh3D.h>

#include <map>

#include <basics/utilities/UbMath.h>

#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbHalfSpace3D.h>

using namespace std;

GbTriangularMesh3D::GbTriangularMesh3D() : GbObject3D()
{
    this->setName("new GbMesh");
    this->nodes     = new vector<GbPoint3D *>;
    this->triangles = new vector<GbTriangle3D *>;
    this->edges     = new vector<GbLine3D *>;

    this->pointinobjecttest = RAYCROSSING;

    this->consistent = false;
    x1min = x1max = x2min = x2max = x3min = x3max = 0.0;
}
/*=============================================================================================*/
GbTriangularMesh3D::GbTriangularMesh3D(string name, vector<GbPoint3D *> *nodes, vector<GbTriangle3D *> *triangles)
    : GbObject3D()
{
    if (name.size() == 0)
        throw UbException(UB_EXARGS, "no name specified");
    if (!nodes)
        throw UbException(UB_EXARGS, "no nodes specified");
    if (!triangles)
        throw UbException(UB_EXARGS, "no triangles specified");

    this->setName(name);
    this->nodes             = nodes;
    this->triangles         = triangles;
    this->edges             = new vector<GbLine3D *>;
    this->pointinobjecttest = RAYCROSSING;

    this->consistent = false;
    x1min = x1max = x2min = x2max = x3min = x3max = 0.0;
}
/*=============================================================================================*/
GbTriangularMesh3D::GbTriangularMesh3D(string name, vector<GbTriangle3D *> *tris) : GbObject3D()
{
    cout << "Das Teil erzeugt seinen KnotenVector aus den Dreiecken ...\n Es sollte deleteRedundantNodes() aufgerufen "
            "werden \n";
    if (name.size() == 0)
        throw UbException(UB_EXARGS, "no name specified");
    if (!tris)
        throw UbException(UB_EXARGS, "no triangles specified");

    vector<GbPoint3D *> *points = new vector<GbPoint3D *>;
    this->triangles             = new vector<GbTriangle3D *>;
    GbPoint3D *p1               = NULL;
    GbPoint3D *p2               = NULL;
    GbPoint3D *p3               = NULL;
    for (int u = 0; u < (int)tris->size(); u++) {
        if (UbMath::zero((*tris)[u]->getArea())) {
            (*tris)[u]->finalize();
            delete (*tris)[u];
            (*tris)[u] = NULL;
            continue;
        }
        this->triangles->push_back((*tris)[u]);
        p1 = (*tris)[u]->getPoint1();
        p2 = (*tris)[u]->getPoint2();
        p3 = (*tris)[u]->getPoint3();
        points->push_back(p1);
        points->push_back(p2);
        points->push_back(p3);
    }

    this->setName(name);
    this->nodes = points;
    // this->triangles        = triangles;
    this->edges = new vector<GbLine3D *>;
    this->edges->resize(0, NULL);
    this->pointinobjecttest = RAYCROSSING;

    this->consistent = false;
    x1min = x1max = x2min = x2max = x3min = x3max = 0.0;
}
/*=============================================================================================*/
GbTriangularMesh3D::GbTriangularMesh3D(string name, vector<GbPoint3D *> *nodes, vector<GbLine3D *> *edges,
                                       vector<GbTriangle3D *> *triangles)
    : GbObject3D()
{
    if (name.size() == 0)
        throw UbException(UB_EXARGS, "no name specified");
    if (!nodes)
        throw UbException(UB_EXARGS, "no nodes specified");
    if (!triangles)
        throw UbException(UB_EXARGS, "no triangles specified");
    if (!edges)
        throw UbException(UB_EXARGS, "no edges specified");

    this->setName(name);
    this->nodes             = nodes;
    this->edges             = edges;
    this->triangles         = triangles;
    this->pointinobjecttest = RAYCROSSING;

    this->consistent = false;
    x1min = x1max = x2min = x2max = x3min = x3max = 0.0;
}
/*=============================================================================================*/
GbTriangularMesh3D::~GbTriangularMesh3D()
{
    if (this->nodes) {
        for (unsigned u = 0; u < nodes->size(); u++)
            delete (*nodes)[u];
        delete nodes;
    }
    if (triangles) {
        for (unsigned u = 0; u < triangles->size(); u++)
            delete (*triangles)[u];
        delete triangles;
    }
}
/*======================================================================*/
void GbTriangularMesh3D::deleteRedundantNodes()
{
    std::map<GbPoint3D *, GbTriangle3D *> pointTriMap;
    GbPoint3D *p1     = NULL;
    GbPoint3D *p2     = NULL;
    GbPoint3D *p3     = NULL;
    GbTriangle3D *tri = NULL;

    for (int u = 0; u < (int)this->triangles->size(); u++) {
        tri = (*this->triangles)[u];
        p1  = tri->getPoint1();
        p2  = tri->getPoint2();
        p3  = tri->getPoint3();
        pointTriMap.insert(pair<GbPoint3D *, GbTriangle3D *>(p1, tri));
        pointTriMap.insert(pair<GbPoint3D *, GbTriangle3D *>(p2, tri));
        pointTriMap.insert(pair<GbPoint3D *, GbTriangle3D *>(p3, tri));
    }

    cout << "Nodes before deleting redundant:" << this->nodes->size() << endl;
    GbPoint3D *pA = NULL;
    GbPoint3D *pB = NULL;
    std::map<GbPoint3D *, GbTriangle3D *>::iterator mapIterator;
    for (int u = 0; u < (int)this->nodes->size(); u++) {
        // cout<<u<<" von "<<this->nodes->size()<<endl;
        pA = (*this->nodes)[u];
        if (pA == NULL)
            continue;
        for (int w = u + 1; w < (int)this->nodes->size(); w++) {
            //   cout<<w<<" Wvon "<<this->nodes->size()<<endl;
            pB = (*this->nodes)[w];
            if (pB == NULL)
                continue;
            if (pA->equals(pB)) {
                // doppelter Knoten ...
                mapIterator = pointTriMap.find(pB);
                tri         = dynamic_cast<GbTriangle3D *>(mapIterator->second);
                if (!tri)
                    throw UbException(UB_EXARGS, "triangle not found");
                p1 = tri->getPoint1();
                p2 = tri->getPoint2();
                p3 = tri->getPoint3();
                if (pB == p1)
                    tri->setPoint(pA, 0);
                else if (pB == p2)
                    tri->setPoint(pA, 1);
                else if (pB == p3)
                    tri->setPoint(pA, 2);
                else
                    throw UbException(UB_EXARGS, "node should be there");
                delete pB;
                (*this->nodes)[w] = NULL;
            }
        }
    }
    vector<GbPoint3D *> *points = new vector<GbPoint3D *>;
    for (int u = 0; u < (int)this->nodes->size(); u++) {
        pA = (*this->nodes)[u];
        if (pA != NULL)
            points->push_back(pA);
    }
    delete this->nodes;
    this->nodes = points;
    cout << "Nodes after deleting redundant:" << this->nodes->size() << endl;

    // nochmal kontrolle ...
    pointTriMap.clear();
    for (int u = 0; u < (int)this->triangles->size(); u++) {
        tri = (*this->triangles)[u];
        p1  = tri->getPoint1();
        p2  = tri->getPoint2();
        p3  = tri->getPoint3();
        pointTriMap.insert(pair<GbPoint3D *, GbTriangle3D *>(p1, tri));
        pointTriMap.insert(pair<GbPoint3D *, GbTriangle3D *>(p2, tri));
        pointTriMap.insert(pair<GbPoint3D *, GbTriangle3D *>(p3, tri));
    }
    for (int u = 0; u < (int)this->nodes->size(); u++) {
        pA = (*this->nodes)[u];
        if (pA == NULL)
            throw UbException(UB_EXARGS, "sollte kein NULL pointer sein ...");
        mapIterator = pointTriMap.find(pA);
        tri         = dynamic_cast<GbTriangle3D *>(mapIterator->second);
        if (!tri)
            throw UbException(UB_EXARGS, "triangle not found");
    }
}
/*======================================================================*/
void GbTriangularMesh3D::translate(const double &x1, const double &x2, const double &x3)
{
    GbPoint3D *pt;
    for (int u = 0; u < (int)this->nodes->size(); u++) {
        pt = (*nodes)[u];
        pt->setX1(pt->getX1Coordinate() + x1);
        pt->setX2(pt->getX2Coordinate() + x2);
        pt->setX3(pt->getX3Coordinate() + x3);
    }
    this->consistent = false;
}
/*======================================================================*/
void GbTriangularMesh3D::rotate(const double &alpha, const double &beta, const double &gamma)
{
    if (!this->consistent)
        this->calculateValues();
    double a1 = this->getX1Centroid();
    double a2 = this->getX2Centroid();
    double a3 = this->getX3Centroid();
    CoordinateTransformation3D trafoFor(a1, a2, a3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    CoordinateTransformation3D trafoBack(a1, a2, a3, 1.0, 1.0, 1.0, alpha, beta, gamma);

    vector<GbPoint3D *> points;
    GbPoint3D *p1 = NULL;
    GbPoint3D *p2 = NULL;
    GbPoint3D *p3 = NULL;
    for (int u = 0; u < (int)this->triangles->size(); u++) {
        p1          = (*triangles)[u]->getPoint1();
        p2          = (*triangles)[u]->getPoint2();
        p3          = (*triangles)[u]->getPoint3();
        double p1x1 = trafoFor.transformForwardToX1Coordinate(p1->x1, p1->x2, p1->x3);
        double p1x2 = trafoFor.transformForwardToX2Coordinate(p1->x1, p1->x2, p1->x3);
        double p1x3 = trafoFor.transformForwardToX3Coordinate(p1->x1, p1->x2, p1->x3);
        double p2x1 = trafoFor.transformForwardToX1Coordinate(p2->x1, p2->x2, p2->x3);
        double p2x2 = trafoFor.transformForwardToX2Coordinate(p2->x1, p2->x2, p2->x3);
        double p2x3 = trafoFor.transformForwardToX3Coordinate(p2->x1, p2->x2, p2->x3);
        double p3x1 = trafoFor.transformForwardToX1Coordinate(p3->x1, p3->x2, p3->x3);
        double p3x2 = trafoFor.transformForwardToX2Coordinate(p3->x1, p3->x2, p3->x3);
        double p3x3 = trafoFor.transformForwardToX3Coordinate(p3->x1, p3->x2, p3->x3);
        p1->x1      = trafoBack.transformBackwardToX1Coordinate(p1x1, p1x2, p1x3);
        p1->x2      = trafoBack.transformBackwardToX2Coordinate(p1x1, p1x2, p1x3);
        p1->x3      = trafoBack.transformBackwardToX3Coordinate(p1x1, p1x2, p1x3);
        p2->x1      = trafoBack.transformBackwardToX1Coordinate(p2x1, p2x2, p2x3);
        p2->x2      = trafoBack.transformBackwardToX2Coordinate(p2x1, p2x2, p2x3);
        p2->x3      = trafoBack.transformBackwardToX3Coordinate(p2x1, p2x2, p2x3);
        p3->x1      = trafoBack.transformBackwardToX1Coordinate(p3x1, p3x2, p3x3);
        p3->x2      = trafoBack.transformBackwardToX2Coordinate(p3x1, p3x2, p3x3);
        p3->x3      = trafoBack.transformBackwardToX3Coordinate(p3x1, p3x2, p3x3);
    }
    this->calculateValues();
}
/*======================================================================*/
/**
 * Returns a string representation of this triangular mesh.
 * @return a string representation of this triangular mesh
 */
string GbTriangularMesh3D::toString()
{
    stringstream ss;
    ss << "GbTriangularMesh3D[";
    ss << (int)this->triangles->size() << "-Triangles, " << (int)this->nodes->size() << "-Nodes, "
       << (int)this->edges->size() << "-Edges" << endl;
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
// string GbTriangularMesh3D::getName(){ return(this->name); }

/**
 * Returns the nodes of this triangular mesh.
 * @return the nodes of this triangular mesh
 */
vector<GbPoint3D *> *GbTriangularMesh3D::getNodes() { return (this->nodes); }
/**
 * Returns the triangles of this triangular mesh.
 * @return the triangles of this triangular mesh
 */
vector<GbTriangle3D *> *GbTriangularMesh3D::getTriangles() { return (this->triangles); }
/**
 * Returns the center x1 coordinate of this triangular mesh.
 * @return the center x1 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX1Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (0.5 * (this->x1min + this->x1max));
}
/**
 * Returns the center x2 coordinate of this triangular mesh.
 * @return the center x2 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX2Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (0.5 * (this->x2min + this->x2max));
}
/**
 * Returns the center x3 coordinate of this triangular mesh.
 * @return the center x3 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX3Centroid()
{
    if (!this->consistent)
        this->calculateValues();
    return (0.5 * (this->x3min + this->x3max));
}

/**
 * Returns the minimum x1 coordinate of this triangular mesh.
 * @return the minimum x1 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX1Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1min);
}
/**
 * Returns the maximum x1 coordinate of this triangular mesh.
 * @return the maximum x1 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX1Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x1max);
}
/**
 * Returns the minimum x2 coordinate of this triangular mesh.
 * @return the minimum x2 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX2Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2min);
}
/**
 * Returns the maximum x2 coordinate of this triangular mesh.
 * @return the maximum x2 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX2Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x2max);
}
/**
 * Returns the minimum x3 coordinate of this triangular mesh.
 * @return the minimum x3 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX3Minimum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3min);
}
/**
 * Returns the maximum x3 coordinate of this triangular mesh.
 * @return the maximum x3 coordinate of this triangular mesh
 */
double GbTriangularMesh3D::getX3Maximum()
{
    if (!this->consistent)
        this->calculateValues();
    return (this->x3max);
}

void GbTriangularMesh3D::calculateValues()
{
    double x1, x2, x3;

    this->x1min = (*this->nodes)[0]->getX1Coordinate();
    this->x1max = (*this->nodes)[0]->getX1Coordinate();
    this->x2min = (*this->nodes)[0]->getX2Coordinate();
    this->x2max = (*this->nodes)[0]->getX2Coordinate();
    this->x3min = (*this->nodes)[0]->getX3Coordinate();
    this->x3max = (*this->nodes)[0]->getX3Coordinate();

    for (int i = 1; i < (int)this->nodes->size(); i++) {
        x1 = (*this->nodes)[i]->getX1Coordinate();
        x2 = (*this->nodes)[i]->getX2Coordinate();
        x3 = (*this->nodes)[i]->getX3Coordinate();
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

/**
 * Returns the total area of this triangular mesh.
 * @return the total area of this triangular mesh
 */
double GbTriangularMesh3D::getArea()
{
    double area = 0.0;
    for (int i = 0; i < (int)this->triangles->size(); i++)
        area += (*this->triangles)[i]->getArea();
    return (area);
}
/**
 * Returns the total volume of this triangular mesh.
 * @return the total volume of this triangular mesh
 */
double GbTriangularMesh3D::getVolume()
{
    GbTriangle3D *triangle;
    GbPoint3D *p1;
    GbPoint3D *p2;
    GbPoint3D *p3;

    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double G3i;
    double volume = 0.0;
    int size      = (int)this->triangles->size();
    for (int u = 0; u < size; u++) {
        triangle = (*this->triangles)[u];
        p1       = triangle->getPoint1();
        p2       = triangle->getPoint2();
        p3       = triangle->getPoint3();
        x1       = p1->getX1Coordinate();
        y1       = p1->getX2Coordinate();
        z1       = p1->getX3Coordinate();
        x2       = p2->getX1Coordinate();
        y2       = p2->getX2Coordinate();
        z2       = p2->getX3Coordinate();
        x3       = p3->getX1Coordinate();
        y3       = p3->getX2Coordinate();
        z3       = p3->getX3Coordinate();
        G3i      = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        volume   = volume + G3i / 6.0;
    }
    return volume;
}
/*===============================================*/
UbTupleDouble3 GbTriangularMesh3D::calculateCenterOfGravity()
{
    GbTriangle3D *triangle;
    GbPoint3D *p1;
    GbPoint3D *p2;
    GbPoint3D *p3;

    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double G3i;
    double rSP1   = 0.0;
    double rSP2   = 0.0;
    double rSP3   = 0.0;
    double volume = 0.0;
    int size      = (int)this->triangles->size();
    for (int u = 0; u < size; u++) {
        triangle = (*this->triangles)[u];
        p1       = triangle->getPoint1();
        p2       = triangle->getPoint2();
        p3       = triangle->getPoint3();
        x1       = p1->getX1Coordinate();
        y1       = p1->getX2Coordinate();
        z1       = p1->getX3Coordinate();
        x2       = p2->getX1Coordinate();
        y2       = p2->getX2Coordinate();
        z2       = p2->getX3Coordinate();
        x3       = p3->getX1Coordinate();
        y3       = p3->getX2Coordinate();
        z3       = p3->getX3Coordinate();
        G3i      = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        volume   = volume + G3i / 6.0;
        rSP1     = rSP1 + G3i * (x1 + x2 + x3);
        rSP2     = rSP2 + G3i * (y1 + y2 + y3);
        rSP3     = rSP3 + G3i * (z1 + z2 + z3);
    }
    rSP1 = rSP1 / (24.0 * volume);
    rSP2 = rSP2 / (24.0 * volume);
    rSP3 = rSP3 / (24.0 * volume);

    return { rSP1, rSP2, rSP3 };
}
/*===============================================*/
UbTupleDouble6 GbTriangularMesh3D::calculateMomentOfInertia(double rhoP)
{
    GbTriangle3D *triangle;
    GbPoint3D *p1;
    GbPoint3D *p2;
    GbPoint3D *p3;

    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double G3i;
    double xx, yy, zz, xy, yz, zx;
    double rSP1   = 0.0;
    double rSP2   = 0.0;
    double rSP3   = 0.0;
    double volume = 0.0;
    double top11  = 0.0;
    double top22  = 0.0;
    double top33  = 0.0;
    double top12  = 0.0;
    double top23  = 0.0;
    double top13  = 0.0;
    int size      = (int)this->triangles->size();
    for (int u = 0; u < size; u++) {
        triangle = (*this->triangles)[u];
        p1       = triangle->getPoint1();
        p2       = triangle->getPoint2();
        p3       = triangle->getPoint3();
        x1       = p1->getX1Coordinate();
        y1       = p1->getX2Coordinate();
        z1       = p1->getX3Coordinate();
        x2       = p2->getX1Coordinate();
        y2       = p2->getX2Coordinate();
        z2       = p2->getX3Coordinate();
        x3       = p3->getX1Coordinate();
        y3       = p3->getX2Coordinate();
        z3       = p3->getX3Coordinate();
        G3i      = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        volume   = volume + G3i / 6.0;
        rSP1     = rSP1 + G3i * (x1 + x2 + x3);
        rSP2     = rSP2 + G3i * (y1 + y2 + y3);
        rSP3     = rSP3 + G3i * (z1 + z2 + z3);
    }
    rSP1 = rSP1 / (24.0 * volume);
    rSP2 = rSP2 / (24.0 * volume);
    rSP3 = rSP3 / (24.0 * volume);

    double x1s = 0.0; // rSP1;//0.0;//
    double x2s = 0.0; // rSP2;//0.0;//
    double x3s = 0.0; // rSP3;//0.0;//

    for (int u = 0; u < size; u++) {
        triangle = (*this->triangles)[u];
        p1       = triangle->getPoint1();
        p2       = triangle->getPoint2();
        p3       = triangle->getPoint3();
        x1       = p1->getX1Coordinate() - x1s;
        y1       = p1->getX2Coordinate() - x2s;
        z1       = p1->getX3Coordinate() - x3s;
        x2       = p2->getX1Coordinate() - x1s;
        y2       = p2->getX2Coordinate() - x2s;
        z2       = p2->getX3Coordinate() - x3s;
        x3       = p3->getX1Coordinate() - x1s;
        y3       = p3->getX2Coordinate() - x2s;
        z3       = p3->getX3Coordinate() - x3s;
        G3i      = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        // rSP1 = rSP1+G3i*(x1+x2+x3)/(24.0*volume);
        // rSP2 = rSP2+G3i*(y1+y2+y3)/(24.0*volume);
        // rSP3 = rSP3+G3i*(z1+z2+z3)/(24.0*volume);
        xx    = x1 * x1 + x2 * x2 + x3 * x3 + x1 * x2 + x2 * x3 + x3 * x1;
        yy    = y1 * y1 + y2 * y2 + y3 * y3 + y1 * y2 + y2 * y3 + y3 * y1;
        zz    = z1 * z1 + z2 * z2 + z3 * z3 + z1 * z2 + z2 * z3 + z3 * z1;
        top11 = top11 + (yy + zz) * rhoP * G3i / 60.;
        top22 = top22 + (xx + zz) * rhoP * G3i / 60.;
        top33 = top33 + (yy + xx) * rhoP * G3i / 60.;
        xy    = 2.0 * (x1 * y1 + x2 * y2 + x3 * y3) + x2 * y3 + x3 * y1 + x1 * y2 + x3 * y2 + x1 * y3 + x2 * y1;
        yz    = 2.0 * (y1 * z1 + y2 * z2 + y3 * z3) + y2 * z3 + y3 * z1 + y1 * z2 + y3 * z2 + y1 * z3 + y2 * z1;
        zx    = 2.0 * (z1 * x1 + z2 * x2 + z3 * x3) + z2 * x3 + z3 * x1 + z1 * x2 + z3 * x2 + z1 * x3 + z2 * x1;
        top12 = top12 - xy * rhoP * G3i / 120.;
        top23 = top23 - yz * rhoP * G3i / 120.;
        top13 = top13 - zx * rhoP * G3i / 120.;
    }
    // Satz von Steiner ...
    top11 = top11 - rhoP * volume * (rSP2 * rSP2 + rSP3 + rSP3);
    top22 = top22 - rhoP * volume * (rSP3 * rSP3 + rSP1 * rSP1);
    top33 = top33 - rhoP * volume * (rSP1 * rSP1 + rSP2 * rSP2);
    top12 = top12 + rhoP * volume * rSP1 * rSP2;
    top23 = top23 + rhoP * volume * rSP2 * rSP3;
    top13 = top13 + rhoP * volume * rSP3 * rSP1;

    cout << "Volume:" << volume << "\n Traegheitsmomente:\n";
    cout << " top11:" << top11 << " top22:" << top22 << " top33:" << top33 << endl;
    cout << " top12:" << top12 << " top23:" << top23 << " top13:" << top13 << endl;

    return { top11, top22, top33, top12, top23, top13 };
}

/**
 * Returns the volume of this triangular mesh within the specified rectangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return the volume of this triangular mesh within the specified rectangle
 * @exception NullPointerException if no triangles are found within the specified rectangle
 */
double GbTriangularMesh3D::getVolumeForRectangle(const double & /*p1x1*/, const double & /*p1x2*/,
                                                 const double & /*p2x1*/, const double & /*p2x2*/)
{
    throw UbException(UB_EXARGS, "not yet implemented");
    //    GbPolygon2D polygon;
    //    double      volume = 0.0;
    //    double      area1 = Math.abs((p1x1-p2x1)*(p1x2-p2x2));
    //    double      area2 = 0.0;
    //    double      t1min, t1max, t2min, t2max;
    //    double      x1, x2;
    //    boolean     f = false;

    //    for(int i=0; i<this.triangles.length; i++)
    //    {
    //    t1min = this.triangles[i].getX1Minimum();
    //    t1max = this.triangles[i].getX1Maximum();
    //    if(GbSystem.less2(t1min, t1max, p1x1, p2x1))    continue;
    //    if(GbSystem.greater2(t1min, t1max, p1x1, p2x1)) continue;

    //    t2min = this.triangles[i].getX2Minimum();
    //    t2max = this.triangles[i].getX2Maximum();
    //    if(GbSystem.less2(t2min, t2max, p1x2, p2x2))    continue;
    //    if(GbSystem.greater2(t2min, t2max, p1x2, p2x2)) continue;

    //    if(GbSystem.inOpenInterval(t1min, p1x1, p2x1) && GbSystem.inOpenInterval(t1max, p1x1, p2x1) &&
    //        GbSystem.inOpenInterval(t2min, p1x2, p2x2) && GbSystem.inOpenInterval(t2max, p1x2, p2x2))
    //    {
    //        volume += this.triangles[i].getVolume();
    //        area2  += this.triangles[i].getArea();
    //        f       = true;
    //    }
    //    else
    //    {
    //        polygon = this.triangles[i].createClippedPolygon3D(p1x1, p1x2, p2x1, p2x2);

    //        if(polygon != null && polygon.size() > 2)
    //        {
    //            try
    //            {
    //                x1      = polygon.getX1Centroid();
    //                x2      = polygon.getX2Centroid();
    //                volume += this.triangles[i].getX3Coordinate(x1, x2) * Math.abs(polygon.getArea());
    //                area2  += Math.abs(polygon.getArea());
    //                f       = true;
    //            }
    //            catch(Exception e){}
    //        }
    //    }
    //    if(GbSystem.greaterEqual(area2, area1)) break;
    //}
    //    if(f) return(volume);
    //    else  throw new NullPointerException();
}

/**
 * Returns the triangles of this triangular mesh located within the specified rectangle (may be an empty array).
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return the triangles of this triangular mesh located within the specified rectangle
 */
vector<GbTriangle3D *> *GbTriangularMesh3D::getTrianglesForRectangle(const double & /*p1x1*/, const double & /*p1x2*/,
                                                                     const double & /*p2x1*/, const double & /*p2x2*/)
{
    throw UbException(UB_EXARGS, "not yet implemented");
    //    QbList      triangleList = new QbList(GbTriangle3D.class);
    //    GbPolygon2D polygon;
    //    double      t1min, t1max, t2min, t2max;
    //    double      area1 = Math.abs((p1x1-p2x1)*(p1x2-p2x2));
    //    double      area2 = 0.0;

    //    for(int i=0; i<this.triangles.length; i++)
    //    {
    //    t1min = this.triangles[i].getX1Minimum();
    //    t1max = this.triangles[i].getX1Maximum();
    //    if(GbSystem.less2(t1min, t1max, p1x1, p2x1))    continue;
    //    if(GbSystem.greater2(t1min, t1max, p1x1, p2x1)) continue;

    //    t2min = this.triangles[i].getX2Minimum();
    //    t2max = this.triangles[i].getX2Maximum();
    //    if(GbSystem.less2(t2min, t2max, p1x2, p2x2))    continue;
    //    if(GbSystem.greater2(t2min, t2max, p1x2, p2x2)) continue;

    //    if(GbSystem.inOpenInterval(t1min, p1x1, p2x1) && GbSystem.inOpenInterval(t1max, p1x1, p2x1) &&
    //        GbSystem.inOpenInterval(t2min, p1x2, p2x2) && GbSystem.inOpenInterval(t2max, p1x2, p2x2))
    //    {
    //        try { triangleList.append(this.triangles[i]); } catch(Exception e){}
    //        area2 += this.triangles[i].getArea();
    //    }
    //    else
    //    {
    //        polygon = this.triangles[i].createClippedPolygon3D(p1x1, p1x2, p2x1, p2x2);
    //        if(polygon != null && polygon.size() > 2)
    //        {
    //            try { triangleList.append(this.triangles[i]); } catch(Exception e){}
    //            area2 += Math.abs(polygon.getArea());
    //        }
    //    }
    //    if(GbSystem.greaterEqual(area2, area1)) break;
    //}
    //    return((GbTriangle3D[])triangleList.getObjectArray());
}
/**
 * Returns the nodes of this triangular mesh located within the specified rectangle (may be an empty array).
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return the nodes of this triangular mesh located within the specified rectangle
 */
vector<GbPoint3D *> *GbTriangularMesh3D::getNodesForRectangle(const double & /*p1x1*/, const double & /*p1x2*/,
                                                              const double & /*p2x1*/, const double & /*p2x2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
    //   QbList nodeList = new QbList(GbPoint3D.class);

    //   for(int i=0; i<this.edges.length; i++)
    //   {
    // if(GbSystem.inClosedInterval(this.nodes[i].getX1Coordinate(), p1x1, p2x1) &&
    // GbSystem.inClosedInterval(this.nodes[i].getX2Coordinate(), p1x2, p2x2))
    //{
    //    try { nodeList.append(this.nodes[i]); } catch(Exception e){}
    //}
    //   }
    //   return((GbPoint3D[])nodeList.getObjectArray());
}

/**
 * Returns the difference of maximum and minimum x3 coordinates
 * of this triangular mesh within the specified rectangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return the difference of maximum and minimum x3 coordinates of this triangular mesh within the specified rectangle
 * @exception NullPointerException if no triangles are found within the specified rectangle
 */
double GbTriangularMesh3D::getX3RangeForRectangle(const double & /*p1x1*/, const double & /*p1x2*/,
                                                  const double & /*p2x1*/, const double & /*p2x2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
    //     GbPolygon3D polygon;
    //     boolean     f     = false;
    //     double      x3min = 0.0;
    //     double      x3max = 0.0;
    //     double      t1min, t1max, t2min, t2max;
    //     double      area1 = Math.abs((p1x1-p2x1)*(p1x2-p2x2));
    //     double      area2 = 0.0;

    //     for(int i=0; i<this.triangles.length; i++)
    //     {
    // t1min = this.triangles[i].getX1Minimum();
    // t1max = this.triangles[i].getX1Maximum();
    // if(GbSystem.less2(t1min, t1max, p1x1, p2x1))    continue;
    // if(GbSystem.greater2(t1min, t1max, p1x1, p2x1)) continue;

    // t2min = this.triangles[i].getX2Minimum();
    // t2max = this.triangles[i].getX2Maximum();
    // if(GbSystem.less2(t2min, t2max, p1x2, p2x2))    continue;
    // if(GbSystem.greater2(t2min, t2max, p1x2, p2x2)) continue;

    // if(GbSystem.inOpenInterval(t1min, p1x1, p2x1) && GbSystem.inOpenInterval(t1max, p1x1, p2x1) &&
    //    GbSystem.inOpenInterval(t2min, p1x2, p2x2) && GbSystem.inOpenInterval(t2max, p1x2, p2x2))
    // {
    //    if(f)
    //    {
    //       if(this.triangles[i].getX3Minimum() < x3min) x3min = this.triangles[i].getX3Minimum();
    //       if(this.triangles[i].getX3Maximum() > x3max) x3max = this.triangles[i].getX3Maximum();
    //    }
    //    else
    //    {
    //       x3min = this.triangles[i].getX3Minimum();
    //       x3max = this.triangles[i].getX3Maximum();
    //       f     = true;
    //    }
    //    area2 += this.triangles[i].getArea();
    //}
    // else
    // {
    //    polygon = this.triangles[i].createClippedPolygon3D(p1x1, p1x2, p2x1, p2x2);

    //    if(polygon != null && polygon.size() > 2)
    //    {
    //       if(f)
    //       {
    //          if(polygon.getX3Minimum() < x3min) x3min = polygon.getX3Minimum();
    //          if(polygon.getX3Maximum() > x3max) x3max = polygon.getX3Maximum();
    //       }
    //       else
    //       {
    //          x3min = polygon.getX3Minimum();
    //          x3max = polygon.getX3Maximum();
    //          f     = true;
    //       }
    //       area2 += Math.abs(polygon.getArea());
    //    }
    // }
    // if(GbSystem.greaterEqual(area2, area1)) break;
    //     }
    //     if(f) return(x3max-x3min);
    //     else  throw new NullPointerException();
}
/**
 * Returns the minimum x3 coordinates of this triangular mesh within the specified rectangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return the minimum x3 coordinates of this triangular mesh within the specified rectangle
 * @exception NullPointerException if no triangles are found within the specified rectangle
 */
double GbTriangularMesh3D::getX3MinimumForRectangle(const double & /*p1x1*/, const double & /*p1x2*/,
                                                    const double & /*p2x1*/, const double & /*p2x2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
    //    GbPolygon3D polygon;
    //    boolean     f     = false;
    //    double      x3min = 0.0;
    //    double      t1min, t1max, t2min, t2max;
    //    double      area1 = Math.abs((p1x1-p2x1)*(p1x2-p2x2));
    //    double      area2 = 0.0;

    //    for(int i=0; i<this.triangles.length; i++)
    //    {
    // t1min = this.triangles[i].getX1Minimum();
    // t1max = this.triangles[i].getX1Maximum();
    // if(GbSystem.less2(t1min, t1max, p1x1, p2x1))    continue;
    // if(GbSystem.greater2(t1min, t1max, p1x1, p2x1)) continue;

    // t2min = this.triangles[i].getX2Minimum();
    // t2max = this.triangles[i].getX2Maximum();
    // if(GbSystem.less2(t2min, t2max, p1x2, p2x2))    continue;
    // if(GbSystem.greater2(t2min, t2max, p1x2, p2x2)) continue;

    // if(GbSystem.inOpenInterval(t1min, p1x1, p2x1) && GbSystem.inOpenInterval(t1max, p1x1, p2x1) &&
    //   GbSystem.inOpenInterval(t2min, p1x2, p2x2) && GbSystem.inOpenInterval(t2max, p1x2, p2x2))
    //{
    //   if(f)
    //   {
    //      if(this.triangles[i].getX3Minimum() < x3min) x3min = this.triangles[i].getX3Minimum();
    //   }
    //   else
    //   {
    //      x3min = this.triangles[i].getX3Minimum();
    //      f     = true;
    //   }
    //   area2 += this.triangles[i].getArea();
    //}
    // else
    //{
    //   polygon = this.triangles[i].createClippedPolygon3D(p1x1, p1x2, p2x1, p2x2);

    //   if(polygon != null && polygon.size() > 2)
    //   {
    //      if(f)
    //      {
    //         if(polygon.getX3Minimum() < x3min) x3min = polygon.getX3Minimum();
    //      }
    //      else
    //      {
    //         x3min = polygon.getX3Minimum();
    //         f     = true;
    //      }
    //      area2 += Math.abs(polygon.getArea());
    //   }
    //}
    // if(GbSystem.greaterEqual(area2, area1)) break;
    //    }
    //    if(f) return(x3min);
    //    else  throw new NullPointerException();
}
/**
 * Returns the maximum x3 coordinates of this triangular mesh within the specified rectangle.
 * @param p1x1 the 1st x1 coordinate of the rectangle
 * @param p1x2 the 1st x2 coordinate of the rectangle
 * @param p2x1 the 2nd x1 coordinate of the rectangle
 * @param p2x2 the 2nd x2 coordinate of the rectangle
 * @return the maximum x3 coordinates of this triangular mesh within the specified rectangle
 * @exception NullPointerException if no triangles are found within the specified rectangle
 */
double GbTriangularMesh3D::getX3MaximumForRectangle(const double & /*p1x1*/, const double & /*p1x2*/,
                                                    const double & /*p2x1*/, const double & /*p2x2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
    //    GbPolygon3D polygon;
    //    boolean     f     = false;
    //    double      x3max = 0.0;
    //    double      t1min, t1max, t2min, t2max;
    //    double      area1 = Math.abs((p1x1-p2x1)*(p1x2-p2x2));
    //    double      area2 = 0.0;

    //    for(int i=0; i<this.triangles.length; i++)
    //    {
    // t1min = this.triangles[i].getX1Minimum();
    // t1max = this.triangles[i].getX1Maximum();
    // if(GbSystem.less2(t1min, t1max, p1x1, p2x1))    continue;
    // if(GbSystem.greater2(t1min, t1max, p1x1, p2x1)) continue;

    // t2min = this.triangles[i].getX2Minimum();
    // t2max = this.triangles[i].getX2Maximum();
    // if(GbSystem.less2(t2min, t2max, p1x2, p2x2))    continue;
    // if(GbSystem.greater2(t2min, t2max, p1x2, p2x2)) continue;

    // if(GbSystem.inOpenInterval(t1min, p1x1, p2x1) && GbSystem.inOpenInterval(t1max, p1x1, p2x1) &&
    //   GbSystem.inOpenInterval(t2min, p1x2, p2x2) && GbSystem.inOpenInterval(t2max, p1x2, p2x2))
    //{
    //   if(f)
    //   {
    //      if(this.triangles[i].getX3Maximum() < x3max) x3max = this.triangles[i].getX3Maximum();
    //   }
    //   else
    //   {
    //      x3max = this.triangles[i].getX3Maximum();
    //      f     = true;
    //   }
    //   area2 += this.triangles[i].getArea();
    //}
    // else
    //{
    //   polygon = this.triangles[i].createClippedPolygon3D(p1x1, p1x2, p2x1, p2x2);

    //   if(polygon != null && polygon.size() > 2)
    //   {
    //      if(f)
    //      {
    //         if(polygon.getX3Maximum() < x3max) x3max = polygon.getX3Maximum();
    //      }
    //      else
    //      {
    //         x3max = polygon.getX3Maximum();
    //         f     = true;
    //      }
    //      area2 += Math.abs(polygon.getArea());
    //   }
    //}
    // if(GbSystem.greaterEqual(area2, area1)) break;
    //    }
    //    if(f) return(x3max);
    //    else  throw new NullPointerException();
}
/*======================================================================*/
vector<GbTriangle3D *> GbTriangularMesh3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> tris;
    GbTriangle3D *triangle;
    GbPoint3D *p1;
    GbPoint3D *p2;
    GbPoint3D *p3;
    int size = (int)this->triangles->size();
    for (int u = 0; u < size; u++) {
        triangle = (*this->triangles)[u];
        p1       = new GbPoint3D(triangle->getPoint1());
        p2       = new GbPoint3D(triangle->getPoint2());
        p3       = new GbPoint3D(triangle->getPoint3());
        tris.push_back(new GbTriangle3D(p1, p2, p3));
    }
    return tris;
}
/*======================================================================*/
/*
 * Function to determine if the point is inside the polyhedron defined as a 3D object
 * using the Halfspace algorithm
 * @param xp the x-coordinate of the point
 * @param yp the y-coordinate of the point
 * @param zp the z-coordinate of the point
 * @return true if point is inside else return false
 */
bool GbTriangularMesh3D::isPointInObject3DHalfSpace(const double &xp, const double &yp, const double &zp)
{
    vector<GbTriangle3D *> *Triangles = this->triangles;
    int Trianglesize                  = (int)Triangles->size();
    // GbPoint3D Point(xp,yp,zp);
    for (int i = 0; i < Trianglesize; i++) {
        // GbPoint3D* point1 = (*Triangles)[i]->getPoint1();
        // GbPoint3D* point2 = (*Triangles)[i]->getPoint2();
        // GbPoint3D* point3 = (*Triangles)[i]->getPoint3();

        // GbHalfSpace3D halfspace(point1, point2, point3);
        GbHalfSpace3D halfspace((*Triangles)[i]);
        if (halfspace.ptInside(xp, yp, zp))
            return false;
    }
    return true;
}
/*======================================================================*/
bool GbTriangularMesh3D::isPointInObject3DRayCrossing(const double &xp, const double &yp, const double &zp, int radius,
                                                      int /*numVertices*/, int numTriangles)
{
    GbVector3D point(xp, yp, zp);

    if (this->InPolyhedron(numTriangles, point, radius))
        return true;
    else
        return false;
}
/*======================================================================*/
bool GbTriangularMesh3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3)
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

    // Halfspace algorithm, Area of spherical polygons algorithm or Ray crossing algorithm
    GbVector3D bMin(boundingCube.getPoint1());
    GbVector3D bMax(boundingCube.getPoint2());

    boundingCube.finalize();

    bMin       = bMax.Subtract(bMin);
    int radius = (int)(bMin.Length() + 1) + 1;

    if (this->pointinobjecttest == HALFSPACE)
        return this->isPointInObject3DHalfSpace(x1, x2, x3);
    else if (this->pointinobjecttest == RAYCROSSING)
        return this->isPointInObject3DRayCrossing(x1, x2, x3, radius, (int)this->nodes->size(),
                                                  (int)this->triangles->size());
    else
        throw UbException(UB_EXARGS, "no ptInObjTest");
}
/*======================================================================*/
bool GbTriangularMesh3D::isPointInGbObject3D(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/,
                                             bool & /*pointIsOnBoundary*/)
{
    throw UbException(UB_EXARGS, "not implemented");
}
/*======================================================================*/
GbLine3D *GbTriangularMesh3D::createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
}
/*======================================================================*/
void GbTriangularMesh3D::writeMesh(string filename, WbWriter *writer, bool writeNormals)
{
    vector<UbTupleFloat3> nodes(triangles->size() * 3);
    vector<UbTupleInt3> tris(triangles->size());

    for (size_t i = 0; i < triangles->size(); i++) {
        GbTriangle3D &tri = *((*triangles)[i]);
        GbPoint3D &node1  = *tri.getPoint(0);
        GbPoint3D &node2  = *tri.getPoint(1);
        GbPoint3D &node3  = *tri.getPoint(2);

        nodes[i * 3] =
            makeUbTuple((float)node1.getX1Coordinate(), (float)node1.getX2Coordinate(), (float)node1.getX3Coordinate());
        nodes[i * 3 + 1] =
            makeUbTuple((float)node2.getX1Coordinate(), (float)node2.getX2Coordinate(), (float)node2.getX3Coordinate());
        nodes[i * 3 + 2] =
            makeUbTuple((float)node3.getX1Coordinate(), (float)node3.getX2Coordinate(), (float)node3.getX3Coordinate());

        tris[i] = makeUbTuple((int)i * 3, (int)i * 3 + 1, (int)i * 3 + 2);
    }
    writer->writeTriangles(filename, nodes, tris);

    if (writeNormals) {
        vector<UbTupleFloat3> lineNodes(triangles->size() * 2);
        vector<UbTupleInt2> lines(triangles->size());
        for (size_t i = 0; i < triangles->size(); i++) {
            GbVector3D vec = (*triangles)[i]->getNormal();
            lineNodes[i * 2] =
                makeUbTuple((float)(*triangles)[i]->getX1Centroid(), (float)(*triangles)[i]->getX2Centroid(),
                            (float)(*triangles)[i]->getX3Centroid());
            lineNodes[i * 2 + 1] = makeUbTuple((float)((*triangles)[i]->getX1Centroid() + vec.X1()),
                                               (float)((*triangles)[i]->getX2Centroid() + vec.X2()),
                                               (float)((*triangles)[i]->getX3Centroid() + vec.X3()));

            lines[i] = makeUbTuple((int)i * 2, (int)i * 2 + 1);
        }
        writer->writeLines(filename + "_normals", lineNodes, lines);
    }
}

/*======================================================================*/
/*
This function returns a char:
'V': the query point a coincides with a Vertex of polyhedron P.
'E': the query point a is in the relative interior of an Edge of polyhedron P.
'F': the query point a is in the relative interior of a Face of polyhedron P.
'i': the query point a is strictly interior to polyhedron P.
'o': the query point a is strictly exterior to( or outside of) polyhedron P.
*/
bool GbTriangularMesh3D::InPolyhedron(int F, GbVector3D &q, int radius)
{
    GbVector3D r; /* Ray endpoint. */
    GbVector3D p; /* Intersection point; not used. */
    int f, k = 0, crossings = 0;
    char code = '?';

    while (k++ < F) {
        crossings = 0;

        RandomRay(r, radius);
        r = q.Add(r);
        // printf("Ray endpoint: (%d,%d,%d)\n", r[0],r[1],r[2] );

        for (f = 0; f < F; f++) /* Begin check each face */
        {
            if (BoxTest((*this->triangles)[f], q, r) == false)
                code = '0'; // printf("BoxTest = 0!\n");
            else
                code = SegTriInt((*this->triangles)[f], q, r,
                                 p); // printf( "Face = %d: BoxTest/SegTriInt returns %c\n\n", f, code );

            /* If ray is degenerate, then goto outer while to generate another. */
            if (code == 'p' || code == 'v' || code == 'e')
                break; // goto LOOP; //printf("Degenerate ray\n");
            /* If ray hits face at interior point, increment crossings. */
            else if (code == 'f')
                crossings++; // printf( "crossings = %d\n", crossings );
            /* If query endpoint q sits on a V/E/F, return that code. */
            else if (code == 'V' || code == 'E' || code == 'F')
                return true;
            /* If ray misses triangle, do nothing. */
            else if (code == '0') { /*nothing to do*/
            } else
                throw UbException(UB_EXARGS, "Error");
        } /* End check each face */

        /* No degeneracies encountered: ray is generic, so finished. */
        if (f >= F)
            break;
    } /* End while loop */

    //   printf( "Crossings = %d\n", crossings );
    /* q strictly interior to polyhedron iff an odd number of crossings. */
    if ((crossings % 2) == 1)
        return true;

    return false;
}

/* Return a random ray endpoint */
void GbTriangularMesh3D::RandomRay(GbVector3D &ray, int radius)
{
    double x, y, z, w, t;

    double MAX_INT = 2147483647;
    /* Generate a random point on a sphere of radius 1. */
    /* the sphere is sliced at z, and a random point at angle t
    generated on the circle of intersection. */
    z = 2.0 * (double)rand() / MAX_INT - 1.0;
    t = 2.0 * UbMath::PI * (double)rand() / MAX_INT;
    w = sqrt(1 - z * z);
    x = w * cos(t);
    y = w * sin(t);

    ray[0] = radius * x;
    ray[1] = radius * y;
    ray[2] = radius * z;

    /*printf( "RandomRay returns %6d %6d %6d\n", ray[X], ray[Y], ray[Z] );*/
}

/*---------------------------------------------------------------------
'p': The segment lies wholly within the plane.
'q': The q endpoint is on the plane (but not 'p').
'r': The r endpoint is on the plane (but not 'p').
'0': The segment lies strictly to one side or the other of the plane.
'1': The segement intersects the plane, and 'p' does not hold.
---------------------------------------------------------------------*/
char GbTriangularMesh3D::SegPlaneInt(GbTriangle3D *T, GbVector3D &q, GbVector3D &r, GbVector3D &p, int *m)
{
    // cout<<"SegPlaneInt..\n";
    GbVector3D N;
    double D;
    GbVector3D rq;
    double num, denom, t;
    int i;

    *m = PlaneCoeff(T, N, &D);
    /*printf("m=%d; plane=(%lf,%lf,%lf,%lf)\n", m, N[X],N[Y],N[Z],D);*/
    num   = D - q.Dot(N);
    rq    = r.Subtract(q);
    denom = rq.Dot(N);
    /*printf("SegPlaneInt: num=%lf, denom=%lf\n", num, denom );*/

    if (denom == 0.0) { /* Segment is parallel to plane. */
        if (num == 0.0) /* q is on plane. */
                        // if (UbMath::zero(denom)) {  /* Segment is parallel to plane. */
            //   if ( UbMath::zero(num))   /* q is on plane. */
            return 'p';
        else
            return '0';
    } else
        t = num / denom;
    /*printf("SegPlaneInt: t=%lf \n", t );*/

    for (i = 0; i < 3; i++)
        p[i] = q[i] + t * (r[i] - q[i]);

    if ((0.0 < t) && (t < 1.0))
        return '1';
    else if (num == 0.0) /* t == 0 */
        return 'q';
    else if (num == denom) /* t == 1 */
        return 'r';
    else
        return '0';

    // if ( (0.0 < t) && (t < 1.0) )
    //   return '1';
    // else if ( UbMath::zero(num))   /* t == 0 */
    //   return 'q';
    // else if ( UbMath::equal(num , denom) ) /* t == 1 */
    //   return 'r';
    // else return '0';
}
/*---------------------------------------------------------------------
Computes N & D and returns index m of largest component.
---------------------------------------------------------------------*/
int GbTriangularMesh3D::PlaneCoeff(GbTriangle3D *T, GbVector3D &N, double *D)
{
    int i;
    double t;             /* Temp storage */
    double biggest = 0.0; /* Largest component of normal vector. */
    int m          = 0;   /* Index of largest component. */

    N = T->getNormal();
    /*printf("PlaneCoeff: N=(%lf,%lf,%lf)\n", N[X],N[Y],N[Z]);*/
    GbVector3D a(T->getPoint1());

    *D = a.Dot(N);

    /* Find the largest component of N. */
    for (i = 0; i < 3; i++) {
        t = std::fabs(N[i]);
        if (t > biggest) {
            biggest = t;
            m       = i;
        }
    }
    return m;
}

/* Assumption: p lies in the plane containing T.
Returns a char:
'V': the query point p coincides with a Vertex of triangle T.
'E': the query point p is in the relative interior of an Edge of triangle T.
'F': the query point p is in the relative interior of a Face of triangle T.
'0': the query point p does not intersect (misses) triangle T.
*/

char GbTriangularMesh3D::InTri3D(GbTriangle3D *T, int m, GbVector3D & /*p*/)
{
    //   int i;           /* Index for X,Y,Z           */
    int j;            /* Index for X,Y             */
                      //   int k;           /* Index for triangle vertex */
    GbVector3D pp;    /* projected p */
    GbVector3D Tp[3]; /* projected T: three new vertices */

    /* Project out coordinate m in both p and the triangular face */
    // j = 0;
    // for ( i = 0; i < 3; i++ ) {
    //   if ( i != m ) {    /* skip largest coordinate */
    //      pp[j] = p[i];
    //      //for ( k = 0; k < 3; k++ )
    //         std::cout<<"aachtung###############################################";
    //      //            Tp[k][j] = Vertices[T[k]][i];
    //      j++;
    //   }
    //}
    j = 0;
    if (m != 0) {
        Tp[0][j] = T->getPoint1()->getX1Coordinate();
        Tp[1][j] = T->getPoint2()->getX1Coordinate();
        Tp[2][j] = T->getPoint3()->getX1Coordinate();
        j++;
    }
    if (m != 1) {
        Tp[0][j] = T->getPoint1()->getX2Coordinate();
        Tp[1][j] = T->getPoint2()->getX2Coordinate();
        Tp[2][j] = T->getPoint3()->getX2Coordinate();
        j++;
    }
    if (m != 2) {
        Tp[0][j] = T->getPoint1()->getX3Coordinate();
        Tp[1][j] = T->getPoint2()->getX3Coordinate();
        Tp[2][j] = T->getPoint3()->getX3Coordinate();
        j++;
    }

    return (InTri2D(Tp, pp));
}

char GbTriangularMesh3D::InTri2D(GbVector3D Tp[3], GbVector3D &pp)
{
    double area0, area1, area2;

    /* compute three AreaSign() values for pp w.r.t. each edge of the face in 2D */
    area0 = AreaSign(pp, Tp[0], Tp[1]);
    area1 = AreaSign(pp, Tp[1], Tp[2]);
    area2 = AreaSign(pp, Tp[2], Tp[0]);
    //  printf("area0=%f  area1=%f  area2=%f\n",area0,area1,area2);

    if (((area0 == 0.) && (area1 > 0.) && (area2 > 0.)) || ((area1 == 0.) && (area0 > 0.) && (area2 > 0.)) ||
        ((area2 == 0.) && (area0 > 0.) && (area1 > 0.)))
        return 'E';

    if (((area0 == 0.) && (area1 < 0.) && (area2 < 0.)) || ((area1 == 0.) && (area0 < 0.) && (area2 < 0.)) ||
        ((area2 == 0.) && (area0 < 0.) && (area1 < 0.)))
        return 'E';

    if (((area0 > 0.) && (area1 > 0.) && (area2 > 0.)) || ((area0 < 0.) && (area1 < 0.) && (area2 < 0.)))
        return 'F';

    if ((area0 == 0.0) && (area1 == 0.0) && (area2 == 0.0))
        fprintf(stderr, "Error in InTriD\n"), exit(EXIT_FAILURE);

    if (((area0 == 0.) && (area1 == 0.)) || ((area0 == 0.) && (area2 == 0.)) || ((area1 == 0.) && (area2 == 0.)))
        return 'V';

    else
        return '0';
}

double GbTriangularMesh3D::AreaSign(GbVector3D &a, GbVector3D &b, GbVector3D &c)
{
    double area2;

    area2 = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);

    return area2;
    /* The area should be an integer. */
    // FIXME: unrechable code
    // if      ( area2 >  0.5 ) return  1;
    // else if ( area2 < -0.5 ) return -1;
    // else                     return  0;
}

char GbTriangularMesh3D::SegTriInt(GbTriangle3D *T, GbVector3D &q, GbVector3D &r, GbVector3D &p)
{
    int code = '?';
    int m    = -1;

    code = (unsigned char)SegPlaneInt(T, q, r, p, &m);
    //  printf("SegPlaneInt code=%c, m=%d; p=(%lf,%lf,%lf)\n", code,m,p[0],p[1],p[2]);

    if (code == '0')
        return '0';
    else if (code == 'q')
        return InTri3D(T, m, q);
    else if (code == 'r')
        return InTri3D(T, m, r);
    else if (code == 'p')
        return InPlane(T, m, q, r, p);
    else if (code == '1')
        return SegTriCross(T, q, r);
    else /* Error */
        return code;
}

char GbTriangularMesh3D::InPlane(GbTriangle3D * /*T*/, int /*m*/, GbVector3D & /*q*/, GbVector3D & /*r*/,
                                 GbVector3D & /*p*/)
{
    // cout<<"inplane\n";
    /* NOT IMPLEMENTED */
    return 'p';
}

/*---------------------------------------------------------------------
The signed volumes of three tetrahedra are computed, determined
by the segment qr, and each edge of the triangle.
Returns a char:
'v': the open segment includes a vertex of T.
'e': the open segment includes a point in the relative interior of an edge
of T.
'f': the open segment includes a point in the relative interior of a face
of T.
'0': the open segment does not intersect triangle T.
---------------------------------------------------------------------*/

char GbTriangularMesh3D::SegTriCross(GbTriangle3D *T, GbVector3D &q, GbVector3D &r)
{
    // cout<<"SegTriCross\n";
    double vol0, vol1, vol2;
    GbVector3D vert0(T->getPoint1());
    GbVector3D vert1(T->getPoint2());
    GbVector3D vert2(T->getPoint3());

    vol0 = VolumeSign(q, vert0, vert1, r);
    vol1 = VolumeSign(q, vert1, vert2, r);
    vol2 = VolumeSign(q, vert2, vert0, r);

    //  printf( "SegTriCross:  vol0 = %d; vol1 = %d; vol2 = %d\n", vol0, vol1, vol2 );

    /* Same sign: segment intersects interior of triangle. */
    if (((vol0 > 0.) && (vol1 > 0.) && (vol2 > 0.)) || ((vol0 < 0.) && (vol1 < 0.) && (vol2 < 0.)))
    // if ( ( UbMath::greater(vol0, 0. ) && UbMath::greater(vol1 , 0. ) && UbMath::greater(vol2 , 0. ) ) ||
    //     ( UbMath::less(vol0, 0. ) && UbMath::less(vol1, 0. ) && UbMath::less(vol2, 0. ) ) )
    {
        return 'f';
    }

    /* Opposite sign: no intersection between segment and triangle */
    if (((vol0 > 0.) || (vol1 > 0.) || (vol2 > 0.)) && ((vol0 < 0.) || (vol1 < 0.) || (vol2 < 0.))) {
        return '0';
    } else if ((vol0 == 0.0) && (vol1 == 0.0) && (vol2 == 0.0)) {
        std::cout << vol0 << " " << vol1 << " " << vol2 << std::endl;
        fprintf(stderr, "Error 1 in SegTriCross\n"), exit(EXIT_FAILURE);
    }

    /* Two zeros: segment intersects vertex. */
    else if (((vol0 == 0.) && (vol1 == 0.)) || ((vol0 == 0.) && (vol2 == 0.)) || ((vol1 == 0.) && (vol2 == 0.))) {
        return 'v';
    }

    /* One zero: segment intersects edge. */
    else if ((vol0 == 0.) || (vol1 == 0.) || (vol2 == 0.)) {
        return 'e';
    }

    throw UbException(UB_EXARGS, "fprintf( stderr, Error 2 in SegTriCross\n ), exit(EXIT_FAILURE);");
}

double GbTriangularMesh3D::VolumeSign(GbVector3D &a, GbVector3D &b, GbVector3D &c, GbVector3D &d)
{
    double vol;
    double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
    double bxdx, bydy, bzdz, cxdx, cydy, czdz;

    ax = a[0];
    ay = a[1];
    az = a[2];
    bx = b[0];
    by = b[1];
    bz = b[2];
    cx = c[0];
    cy = c[1];
    cz = c[2];
    dx = d[0];
    dy = d[1];
    dz = d[2];

    bxdx = bx - dx;
    bydy = by - dy;
    bzdz = bz - dz;
    cxdx = cx - dx;
    cydy = cy - dy;
    czdz = cz - dz;
    vol  = (az - dz) * (bxdx * cydy - bydy * cxdx) + (ay - dy) * (bzdz * cxdx - bxdx * czdz) +
          (ax - dx) * (bydy * czdz - bzdz * cydy);

    // std::cout<< vol<<std::endl;
    return vol;

    /* The volume should be an integer. */
    // FIXME: unreachable code
    // if      ( vol > 0.5 )   return  1;
    // else if ( vol < -0.5 )  return -1;
    // else                    return  0;
}

bool GbTriangularMesh3D::BoxTest(GbTriangle3D *triangle, GbVector3D &PointQ, GbVector3D &PointR)
{
    double minX1 = triangle->getX1Minimum();
    double minX2 = triangle->getX2Minimum();
    double minX3 = triangle->getX3Minimum();

    double maxX1 = triangle->getX1Maximum();
    double maxX2 = triangle->getX2Maximum();
    double maxX3 = triangle->getX3Maximum();

    if ((PointQ.X1() < minX1) && (PointR.X1() < minX1))
        return false;
    if ((PointQ.X2() < minX2) && (PointR.X2() < minX2))
        return false;
    if ((PointQ.X3() < minX3) && (PointR.X3() < minX3))
        return false;
    if ((PointQ.X1() > maxX1) && (PointR.X1() > maxX1))
        return false;
    if ((PointQ.X2() > maxX2) && (PointR.X2() > maxX2))
        return false;
    if ((PointQ.X3() > maxX3) && (PointR.X3() > maxX3))
        return false;

    return true;
}

//! \}
