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
#include <geometry3d/GbTriFaceMesh3D.h>

#include <basics/Timer/Timer.h>
#include <basics/utilities/UbFileInputASCII.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbRandom.h>
#include <basics/writer/WbWriter.h>

#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbHalfSpace3D.h>

#include <geometry3d/KdTree/KdTree.h>
#include <geometry3d/KdTree/intersectionhandler/KdCountLineIntersectionHandler.h>
#include <geometry3d/KdTree/intersectionhandler/KdCountRayIntersectionHandler.h>
#include <geometry3d/KdTree/splitalgorithms/KdSAHSplit.h>
#include <geometry3d/KdTree/splitalgorithms/KdSpatiallMedianSplit.h>

#define MAX_ITER 10

using namespace std;

GbTriFaceMesh3D::GbTriFaceMesh3D() : GbObject3D()
{
    this->setName("CAB_GbTriFaceMesh3D");
    this->nodes          = new vector<Vertex>;
    this->triangles      = new vector<TriFace>;
    this->consistent     = false;
    this->kdtreeSplitAlg = KDTREE_SAHPLIT;
}
/*=======================================================================*/
GbTriFaceMesh3D::GbTriFaceMesh3D(string name, vector<Vertex> *nodes, vector<TriFace> *triangles,
                                 KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
    : GbObject3D(), nodes(nodes), triangles(triangles), kdtreeSplitAlg(splitAlg)
{
    if (name.empty())
        throw UbException(UB_EXARGS, "no name specified");
    if (!nodes)
        throw UbException(UB_EXARGS, "no nodes specified");
    if (!triangles)
        throw UbException(UB_EXARGS, "no triangles specified");

    this->setName(name);

    if (removeRedundantNodes) {
        this->deleteRedundantNodes(); // dort wird autoamtisch calculateValues() aufgerufen
    } else {
        this->calculateValues();
    }
}
/*=======================================================================*/
GbTriFaceMesh3D::~GbTriFaceMesh3D()
{
    if (nodes) {
        delete nodes;
        nodes = NULL;
    }
    if (triangles) {
        delete triangles;
        triangles = NULL;
    }
    if (kdTree) {
        delete kdTree;
        kdTree = NULL;
    }
}
/*======================================================================*/
void GbTriFaceMesh3D::init()
{
    //nodes      = NULL;
    //triangles  = NULL;
    x1min      = 0.0;
    x1max      = 0.0;
    x1center   = 0.0;
    x2min      = 0.0;
    x2max      = 0.0;
    x2center   = 0.0;
    x3min      = 0.0;
    x3max      = 0.0;
    x3center   = 0.0;
    consistent = false;
}
/*======================================================================*/
GbTriFaceMesh3D *GbTriFaceMesh3D::clone()
{
    vector<GbTriFaceMesh3D::Vertex> *newNodes      = new vector<GbTriFaceMesh3D::Vertex>;
    vector<GbTriFaceMesh3D::TriFace> *newTriangles = new vector<GbTriFaceMesh3D::TriFace>;

    int numberNodes = (int)this->nodes->size();

    double x, y, z;
    for (int u = 0; u < numberNodes; u++) {
        x = (*nodes)[u].x;
        y = (*nodes)[u].y;
        z = (*nodes)[u].z;
        newNodes->push_back(GbTriFaceMesh3D::Vertex((float)x, (float)y, (float)z));
    }
    int numberTris = (int)this->triangles->size();
    UBLOG(logDEBUG1, "numberTris:" << numberTris);

    int id1, id2, id3;
    for (int u = 0; u < numberTris; u++) {
        id1 = (*this->triangles)[u].v1;
        id2 = (*this->triangles)[u].v2;
        id3 = (*this->triangles)[u].v3;
        newTriangles->push_back(GbTriFaceMesh3D::TriFace(id1, id2, id3));
        // cout<<u<<" - id1,id2,id3:"<<id1<<","<<id2<<","<<id3<<endl;
    }
    UBLOG(logDEBUG1, "Tris gelesen");

    GbTriFaceMesh3D *mesh = new GbTriFaceMesh3D("no name", newNodes, newTriangles);
    UBLOG(logDEBUG1, "mesh cloned ...");

    return mesh;
}

/*======================================================================*/
// checks for doppelt nodes und fixed Dreicke die zweimal denselben Knoten haben
void GbTriFaceMesh3D::deleteRedundantNodes()
{
    UBLOG(logDEBUG1,
          "GbTriFaceMesh3D::deleteRedundantNodes - Nodes before deleting redundant: " << this->nodes->size());

    map<Vertex, size_t /*new vecIndex*/> vertexMap;
    map<Vertex, size_t /*new vecIndex*/>::iterator pos;
    map<Vertex, size_t /*new vecIndex*/>::iterator it;

    vector<TriFace> &tris    = *this->triangles;
    vector<Vertex> &oldNodes = *this->nodes;
    vector<Vertex> newNodes;

    for (size_t t = 0; t < tris.size(); t++) {
        if (t % 100 == 0) {
            UBLOG(logDEBUG5, "GbTriFaceMesh3D::deleteRedundantNodes - tri: " << (t) << " von " << tris.size());
        }
        TriFace &tri = tris[t];
        // Knoten bereits in neuem node vector?
        for (int v = 0; v <= 2; v++) {
            Vertex &vert = tri.getNode(v, oldNodes);
            // pos=vertexMap.find( vert );
            // if( pos==vertexMap.end() )
            {
                for (pos = vertexMap.begin(); pos != vertexMap.end(); pos++) {
                    Vertex rhs = pos->first;
                    // if(UbMath::inClosedInterval(vert.z,0.01999, 0.02001))
                    if (fabs(vert.x - rhs.x) < 1.E-5 && fabs(vert.y - rhs.y) < 1.E-5 && fabs(vert.z - rhs.z) < 1.E-5) {
                        break;
                    }
                }
            }
            if (pos != vertexMap.end())
                tri.setNode(v, (int)pos->second);
            else {
                newNodes.push_back(vert);
                int index       = (int)newNodes.size() - 1;
                vertexMap[vert] = index;
                tri.setNode(v, index);
            }
        }
    }

    std::swap(*nodes, newNodes);

    UBLOG(logDEBUG1, "GbTriFaceMesh3D::deleteRedundantNodes - Nodes after deleting redundant:" << this->nodes->size());
    //
    // Das geht irgendwie nicht ...
    //
    // UBLOG(logDEBUG1,"GbTriFaceMesh3D::deleteRedundantNodes - checking for double triangles !!!");
    // UBLOG(logDEBUG1,"GbTriFaceMesh3D::deleteRedundantNodes - Triangles before deleting redundant:
    // "<<this->triangles->size()); vector<TriFace> newSingleTris; newSingleTris.reserve( this->triangles->size() );
    // for(size_t t=0; t<tris.size(); t++)
    //{
    //   Vertex& v1 = tris[t].getNode(0,*nodes);
    //   Vertex& v2 = tris[t].getNode(1,*nodes);
    //   Vertex& v3 = tris[t].getNode(2,*nodes);

    //   if(UbMath::greater(std::fabs(v1.x), 0.0634) && UbMath::inClosedInterval(v1.z, 0.01999, 0.02001))
    //   {
    //      UBLOG2(logINFO,std::cout, "V1:"<<v1.x<<" "<<v1.y<<" "<<v1.z);
    //   }
    //   if(UbMath::greater(std::fabs(v2.x), 0.0634) && UbMath::inClosedInterval(v2.z, 0.01999, 0.02001))
    //   {
    //      UBLOG2(logINFO,std::cout, "V2:"<<v2.x<<" "<<v2.y<<" "<<v2.z);
    //   }
    //   if(UbMath::greater(std::fabs(v3.x), 0.0634) && UbMath::inClosedInterval(v3.z, 0.01999, 0.02001))
    //   {
    //      UBLOG2(logINFO,std::cout, "V3:"<<v3.x<<" "<<v3.y<<" "<<v3.z);
    //   }

    //   bool inList = false;
    //   for(size_t u=0; u<newSingleTris.size(); u++)
    //   {
    //      Vertex& vn1 = newSingleTris[t].getNode(0,*nodes);
    //      Vertex& vn2 = newSingleTris[t].getNode(1,*nodes);
    //      Vertex& vn3 = newSingleTris[t].getNode(2,*nodes);

    //      if(v1==vn1 && v2==vn2 && v3==vn3)      inList = true;
    //      else if(v1==vn1 && v2==vn3 && v3==vn2) inList = true;
    //      else if(v1==vn2 && v2==vn3 && v3==vn1) inList = true;
    //      else if(v1==vn2 && v2==vn1 && v3==vn3) inList = true;
    //      else if(v1==vn3 && v2==vn1 && v3==vn2) inList = true;
    //      else if(v1==vn3 && v2==vn2 && v3==vn1) inList = true;
    //   }
    //   if(!inList) newSingleTris.push_back(tris[t]);
    //   else
    //      UBLOG(logDEBUG1,"GbTriFaceMesh3D::deleteRedundantNodes - inList !!!!");

    //}
    // swap(tris,newSingleTris);

    // UBLOG(logDEBUG1,"GbTriFaceMesh3D::deleteRedundantNodes - Triangles after deleting
    // redundant:"<<this->triangles->size());
    UBLOG(logDEBUG1, "GbTriFaceMesh3D::deleteRedundantNodes - checking for triangles that have same node several times "
                     "or are lines!!!");
    int counter1 = 0;
    int counter2 = 0;
    vector<TriFace> newTris;
    newTris.reserve(this->triangles->size());
    for (size_t t = 0; t < tris.size(); t++) {
        Vertex &v1 = tris[t].getNode(0, *nodes);
        Vertex &v2 = tris[t].getNode(1, *nodes);
        Vertex &v3 = tris[t].getNode(2, *nodes);
        if (v1 == v2 || v1 == v3 || v2 == v3) {
            counter1++;
        } else if (tris[t].getArea(*nodes) < 1.0E-8) {
            counter2++;
        } else
            newTris.push_back(tris[t]);
    }
    if (counter1) {
        UBLOG(logDEBUG1, "GbTriFaceMesh3D::deleteRedundantNodes - ### Warning ###: found and removed  "
                             << counter1 << " triangle with double nodes!");
    }
    if (counter2) {
        UBLOG(logDEBUG1, "GbTriFaceMesh3D::deleteRedundantNodes - ### Warning ###: found and removed  "
                             << counter2 << " triangle that are lines!");
    }
    if (!counter1 && !counter2) {
        UBLOG(logDEBUG1, "GbTriFaceMesh3D::deleteRedundantNodes - alles gut... nix doppelt");
    } else
        swap(tris, newTris);

    UBLOG(logDEBUG1, "GbTriFaceMesh3D::deleteRedundantNodes - done");
    this->calculateValues();
}
/*======================================================================*/
void GbTriFaceMesh3D::setKdTreeSplitAlgorithm(KDTREE_SPLITAGORITHM mode)
{
    if (kdTree && mode != this->kdtreeSplitAlg) {
        delete kdTree;
        kdTree = NULL;
    }
    this->kdtreeSplitAlg = mode;
}
/*======================================================================*/
/**
 * Returns a string representation of this triangular mesh.
 * @return a string representation of this triangular mesh
 */
string GbTriFaceMesh3D::toString()
{
    stringstream ss;
    ss << "GbTriFaceMesh3D[";
    ss << (int)this->triangles->size() << "-Triangles, " << (int)this->nodes->size() << "-Nodes, " << endl;
    ss << "]";
    return (ss.str());
}
/**
 * Returns the nodes of this triangular mesh.
 * @return the nodes of this triangular mesh
 */
vector<GbTriFaceMesh3D::Vertex> *GbTriFaceMesh3D::getNodes() { return this->nodes; }
/**
 * Returns the triangles of this triangular mesh.
 * @return the triangles of this triangular mesh
 */
vector<GbTriFaceMesh3D::TriFace> *GbTriFaceMesh3D::getTriangles() { return this->triangles; }
/**
 * Returns the center x1 coordinate of this triangular mesh.
 * @return the center x1 coordinate of this triangular mesh
 */
double GbTriFaceMesh3D::getVolume()
{
    vector<Vertex> &vertices = *nodes;
    vector<TriFace> &tris    = *triangles;

    double x1, x2, x3, y1, y2, y3, z1, z2, z3, G3i;
    // double rSP1 = 0.0;double rSP2 = 0.0;double rSP3 = 0.0;
    double volume = 0.0;
    for (size_t t = 0; t < tris.size(); t++) {
        TriFace &triangle = tris[t];
        x1                = triangle.getV1x(vertices);
        y1                = triangle.getV1y(vertices);
        z1                = triangle.getV1z(vertices);
        x2                = triangle.getV2x(vertices);
        y2                = triangle.getV2y(vertices);
        z2                = triangle.getV2z(vertices);
        x3                = triangle.getV3x(vertices);
        y3                = triangle.getV3y(vertices);
        z3                = triangle.getV3z(vertices);
        G3i               = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        volume            = volume + G3i / 6.0;
    }
    return volume;
}
/*===============================================*/
UbTupleDouble3 GbTriFaceMesh3D::calculateCenterOfGravity()
{
    vector<Vertex> &vertices = *nodes;
    vector<TriFace> &tris    = *triangles;

    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double G3i;
    double rSP1 = 0.0, rSP2 = 0.0, rSP3 = 0.0, volume = 0.0;

    for (size_t t = 0; t < tris.size(); t++) {
        TriFace &triangle = tris[t];
        x1                = triangle.getV1x(vertices);
        y1                = triangle.getV1y(vertices);
        z1                = triangle.getV1z(vertices);
        x2                = triangle.getV2x(vertices);
        y2                = triangle.getV2y(vertices);
        z2                = triangle.getV2z(vertices);
        x3                = triangle.getV3x(vertices);
        y3                = triangle.getV3y(vertices);
        z3                = triangle.getV3z(vertices);
        G3i               = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        volume            = volume + G3i / 6.0;
        rSP1              = rSP1 + G3i * (x1 + x2 + x3);
        rSP2              = rSP2 + G3i * (y1 + y2 + y3);
        rSP3              = rSP3 + G3i * (z1 + z2 + z3);
    }
    rSP1 = rSP1 / (24.0 * volume);
    rSP2 = rSP2 / (24.0 * volume);
    rSP3 = rSP3 / (24.0 * volume);

    return { rSP1, rSP2, rSP3 };
}
/*===============================================*/
UbTupleDouble6 GbTriFaceMesh3D::calculateMomentOfInertia(double rhoP)
{
    vector<Vertex> &vertices = *nodes;

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
        TriFace &triangle = (*this->triangles)[u];
        x1                = triangle.getV1x(vertices);
        y1                = triangle.getV1y(vertices);
        z1                = triangle.getV1z(vertices);
        x2                = triangle.getV2x(vertices);
        y2                = triangle.getV2y(vertices);
        z2                = triangle.getV2z(vertices);
        x3                = triangle.getV3x(vertices);
        y3                = triangle.getV3y(vertices);
        z3                = triangle.getV3z(vertices);
        G3i               = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
        volume            = volume + G3i / 6.0;
        rSP1              = rSP1 + G3i * (x1 + x2 + x3);
        rSP2              = rSP2 + G3i * (y1 + y2 + y3);
        rSP3              = rSP3 + G3i * (z1 + z2 + z3);
    }
    rSP1 = rSP1 / (24.0 * volume);
    rSP2 = rSP2 / (24.0 * volume);
    rSP3 = rSP3 / (24.0 * volume);

    double x1s = 0.0; // rSP1;//0.0;//
    double x2s = 0.0; // rSP2;//0.0;//
    double x3s = 0.0; // rSP3;//0.0;//

    for (int u = 0; u < size; u++) {
        TriFace &triangle = (*this->triangles)[u];
        x1                = triangle.getV1x(vertices) - x1s;
        y1                = triangle.getV1y(vertices) - x2s;
        z1                = triangle.getV1z(vertices) - x3s;
        x2                = triangle.getV2x(vertices) - x1s;
        y2                = triangle.getV2y(vertices) - x2s;
        z2                = triangle.getV2z(vertices) - x3s;
        x3                = triangle.getV3x(vertices) - x1s;
        y3                = triangle.getV3y(vertices) - x2s;
        z3                = triangle.getV3z(vertices) - x3s;
        G3i               = x1 * (y2 * z3 - z2 * y3) + y1 * (z2 * x3 - x2 * z3) + z1 * (x2 * y3 - y2 * x3);
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
/*==============================================================*/
void GbTriFaceMesh3D::calculateValues()
{
    relationVertTris.clear();

    if (nodes->empty()) {
        x1min = x1max = x2min = x2max = x3min = x3max = 0.0;
    } else {
        Vertex &v = (*nodes)[0];
        x1min = x1max = v.x;
        x2min = x2max = v.y;
        x3min = x3max = v.z;

        for (size_t i = 1; i < this->nodes->size(); i++) {
            Vertex &v1 = (*nodes)[i];

            x1min = UbMath::min<double>(x1min, v1.x);
            x2min = UbMath::min<double>(x2min, v1.y);
            x3min = UbMath::min<double>(x3min, v1.z);

            x1max = UbMath::max<double>(x1max, v1.x);
            x2max = UbMath::max<double>(x2max, v1.y);
            x3max = UbMath::max<double>(x3max, v1.z);
        }
        x1center = 0.5 * (x1min + x1max);
        x2center = 0.5 * (x2min + x2max);
        x3center = 0.5 * (x3min + x3max);

        vector<TriFace> &tris = *this->triangles;
        vector<Vertex> &verts = *this->nodes;
        for (size_t i = 0; i < this->triangles->size(); i++) {
            tris[i].calculateNormal(verts);
        }
        // relation Vertex <-> Triangle ermitteln
        if (buildVertTriRelationMap) {
            for (size_t t = 0; t < tris.size(); t++) {
                TriFace &tri = tris[t];
                relationVertTris.insert(make_pair(&verts[tri.v1], &tri));
                relationVertTris.insert(make_pair(&verts[tri.v2], &tri));
                relationVertTris.insert(make_pair(&verts[tri.v3], &tri));
            }
        }
    }
    if (kdTree) {
        delete kdTree;
        kdTree = NULL;
    }

    this->consistent = true;
}
/*=========================================================================*/
std::vector<GbTriFaceMesh3D::TriFace *> GbTriFaceMesh3D::getTrianglesForVertex(Vertex *vertex)
{
    if (!buildVertTriRelationMap) {
        buildVertTriRelationMap = true;
        consistent              = false;
    }
    if (!consistent)
        this->calculateValues();

    typedef std::multimap<Vertex *, TriFace *>::iterator Iterator;
    pair<Iterator, Iterator> objRange = relationVertTris.equal_range(vertex);

    std::vector<TriFace *> tmpTris;
    for (Iterator pos = objRange.first; pos != objRange.second; ++pos)
        tmpTris.push_back(pos->second);

    return tmpTris;
}
/*=======================================================*/
void GbTriFaceMesh3D::setCenterCoordinates(const double &x1, const double &x2, const double &x3)
{
    this->translate(x1 - getX1Centroid(), x2 - getX2Centroid(), x3 - getX3Centroid());
}

/*======================================================================*/
void GbTriFaceMesh3D::setCenterCoordinates(const UbTupleDouble3 & /*position*/)
{
    throw UbException(UB_EXARGS, "not implemented for " + (std::string) typeid(*this).name());
}

/*======================================================================*/
void GbTriFaceMesh3D::scale(const double &sx1, const double &sx2, const double &sx3)
{
    CoordinateTransformation3D trafoForw(this->getX1Centroid(), this->getX2Centroid(), this->getX3Centroid(), 1.0, 1.0,
                                         1.0, 0.0, 0.0, 0.0);
    CoordinateTransformation3D trafoBack(this->getX1Centroid(), this->getX2Centroid(), this->getX3Centroid(), sx1, sx2,
                                         sx3, 0, 0, 0);

    vector<Vertex> &vertices = *nodes;
    for (size_t i = 0; i < vertices.size(); i++) {
        Vertex &v   = vertices[i];
        double p1x1 = trafoForw.transformForwardToX1Coordinate(v.x, v.y, v.z);
        double p1x2 = trafoForw.transformForwardToX2Coordinate(v.x, v.y, v.z);
        double p1x3 = trafoForw.transformForwardToX3Coordinate(v.x, v.y, v.z);
        v.x         = (float)trafoBack.transformBackwardToX1Coordinate(p1x1, p1x2, p1x3);
        v.y         = (float)trafoBack.transformBackwardToX2Coordinate(p1x1, p1x2, p1x3);
        v.z         = (float)trafoBack.transformBackwardToX3Coordinate(p1x1, p1x2, p1x3);
    }
    this->calculateValues();
}
/*======================================================================*/
void GbTriFaceMesh3D::rotate(const double &alpha, const double &beta, const double &gamma)
{
    CoordinateTransformation3D trafoForw(this->getX1Centroid(), this->getX2Centroid(), this->getX3Centroid(), 1.0, 1.0,
                                         1.0, 0.0, 0.0, 0.0);
    CoordinateTransformation3D trafoBack(this->getX1Centroid(), this->getX2Centroid(), this->getX3Centroid(), 1.0, 1.0,
                                         1.0, alpha, beta, gamma);

    vector<Vertex> &vertices = *nodes;
    for (size_t i = 0; i < vertices.size(); i++) {
        Vertex &v   = vertices[i];
        double p1x1 = trafoForw.transformForwardToX1Coordinate(v.x, v.y, v.z);
        double p1x2 = trafoForw.transformForwardToX2Coordinate(v.x, v.y, v.z);
        double p1x3 = trafoForw.transformForwardToX3Coordinate(v.x, v.y, v.z);
        v.x         = (float)trafoBack.transformBackwardToX1Coordinate(p1x1, p1x2, p1x3);
        v.y         = (float)trafoBack.transformBackwardToX2Coordinate(p1x1, p1x2, p1x3);
        v.z         = (float)trafoBack.transformBackwardToX3Coordinate(p1x1, p1x2, p1x3);
    }
    this->calculateValues();
}
/*======================================================================*/
void GbTriFaceMesh3D::rotateAroundPoint(const double &px1, const double &px2, const double &px3, const double &alpha,
                                        const double &beta, const double &gamma)
{
    CoordinateTransformation3D trafoForw(px1, px2, px3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    CoordinateTransformation3D trafoBack(px1, px2, px3, 1.0, 1.0, 1.0, alpha, beta, gamma);

    vector<Vertex> &vertices = *nodes;
    for (size_t i = 0; i < vertices.size(); i++) {
        Vertex &v   = vertices[i];
        double p1x1 = trafoForw.transformForwardToX1Coordinate(v.x, v.y, v.z);
        double p1x2 = trafoForw.transformForwardToX2Coordinate(v.x, v.y, v.z);
        double p1x3 = trafoForw.transformForwardToX3Coordinate(v.x, v.y, v.z);
        v.x         = (float)trafoBack.transformBackwardToX1Coordinate(p1x1, p1x2, p1x3);
        v.y         = (float)trafoBack.transformBackwardToX2Coordinate(p1x1, p1x2, p1x3);
        v.z         = (float)trafoBack.transformBackwardToX3Coordinate(p1x1, p1x2, p1x3);
    }
    this->calculateValues();
}

/*======================================================================*/
void GbTriFaceMesh3D::translate(const double &x1, const double &x2, const double &x3)
{
    vector<Vertex> &vertices = *nodes;
    for (size_t i = 0; i < vertices.size(); i++) {
        Vertex &v = vertices[i];
        v.x += static_cast<float>(x1);
        v.y += static_cast<float>(x2);
        v.z += static_cast<float>(x3);
    }
    this->calculateValues();
}
/*======================================================================*/
vector<GbTriangle3D *> GbTriFaceMesh3D::getSurfaceTriangleSet()
{
    // SirAnn: eine miese Speicherlochmethode
    //        hier werden dynmamische Objekte angelegt
    //        mit sowas rechnet von aussen kein Mensch!!!
    vector<GbTriangle3D *> tris(triangles->size());

    for (size_t i = 0; i < this->triangles->size(); i++) {
        Vertex &v1 = (*nodes)[(*triangles)[i].v1];
        Vertex &v2 = (*nodes)[(*triangles)[i].v2];
        Vertex &v3 = (*nodes)[(*triangles)[i].v3];

        tris[i] = new GbTriangle3D(new GbPoint3D(v1.x, v1.y, v1.z), new GbPoint3D(v2.x, v2.y, v2.z),
                                   new GbPoint3D(v3.x, v3.y, v3.z));
    }
    return tris;
}
/*=======================================================*/
void GbTriFaceMesh3D::addSurfaceTriangleSet(vector<UbTupleFloat3> &pts, vector<UbTupleInt3> &tris)
{
    int nodeNr = int(pts.size());
    for (int i = 0; i < (int)this->triangles->size(); i++) {
        Vertex &v1 = (*nodes)[(*triangles)[i].v1];
        Vertex &v2 = (*nodes)[(*triangles)[i].v2];
        Vertex &v3 = (*nodes)[(*triangles)[i].v3];
        pts.push_back(makeUbTuple(v1.x, v1.y, v1.z));
        pts.push_back(makeUbTuple(v2.x, v2.y, v2.z));
        pts.push_back(makeUbTuple(v3.x, v3.y, v3.z));

        tris.push_back(makeUbTuple(nodeNr, nodeNr + 1, nodeNr + 2));
        nodeNr += 3;
    }
}
/*======================================================================*/
// bool GbTriFaceMesh3D::isPointInGbObject3D(const double& x1, const double& x2, const double& x3, int counter)
//{
//
//
//   if( !nodes->empty() )
//   {
//      //Baum erstellen, wen noch keiner vorhanden
//      if( !kdTree)
//      {
//         UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree start");
//         vf::basics::Timer timer; timer.start();
//         if(kdtreeSplitAlg == KDTREE_SAHPLIT     )
//         {
//            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit");
//            this->kdTree = new Kd::Tree<double>( *this, Kd::SAHSplit<double>()            );
//         }
//         else if(kdtreeSplitAlg == KDTREE_SPATIALSPLIT)
//         {
//            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SpatialMedianSplit");
//            this->kdTree = new Kd::Tree<double>( *this, Kd::SpatialMedianSplit<double>() );
//         }
//         else throw UbException(UB_EXARGS, "unknown kdtree split option)" );
//         UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - built kdTree in "<<timer.getCurrentRuntimeInSeconds()<<"seconds");
//      }
//
//      //eigentlicher PIO-Test
//      //int iSec;
//      //for(int i=0; i<100; i++)
//      //{
//      //   Kd::Ray<double> ray(  x1, x2, x3  //, 1, 0 ,0 );
//      //                        , ( x1 < x1center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
//      //                        , ( x2 < x2center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
//      //                        , ( x3 < x3center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
//      );
//      //
//      //   iSec = kdTree->intersectRay( ray, Kd::CountRayIntersectionHandler<double>() );
//      //
//      //   if( iSec != Kd::Intersection::INTERSECT_EDGE ) //KEINE Kante getroffen
//      //   {
//      //      if(iSec == Kd::Intersection::ON_BOUNDARY )
//      //      {
//      //         return true;
//      //      }
//      //      return (iSec&1);  //ungerade anzahl an schnitten --> drinnen
//      //   }
//      //   UBLOG(logDEBUG3, "GbTriFaceMesh3D.isPointInGbObject3D.if  - an edge was hit ");
//      //}
//      //throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
//      int iSec1,iSec2;
//
//      Kd::Ray<double> ray1(  x1, x2, x3, 1.0, 0.0 ,0.0 );
//      iSec1 = kdTree->intersectRay( ray1, Kd::CountRayIntersectionHandler<double>() );
//      Kd::Ray<double> ray2(  x1, x2, x3, -1.0, 0.0 ,0.0 );
//      iSec2 = kdTree->intersectRay( ray2, Kd::CountRayIntersectionHandler<double>() );
//
//      if(iSec1 == Kd::Intersection::ON_BOUNDARY || iSec2 == Kd::Intersection::ON_BOUNDARY)
//      {
//         return true;
//      }
//      if( iSec1 == Kd::Intersection::INTERSECT_EDGE && iSec2 == Kd::Intersection::INTERSECT_EDGE)
//      {
//         UBLOG(logINFO, "GbTriFaceMesh3D.isPointInGbObject3D.INTERSECT_EDGE");
//         double eps = UbMath::getEqualityEpsilon<float>()*1000.0;
//         if (counter>100) {return(iSec1&1);  UBLOG(logINFO, "NACH 100 Iterationen Eps umsetzen aufgegeben!");}
//         return this->isPointInGbObject3D(x1+eps, x2+eps, x3+eps,(counter+1));
//      }
//      else if( iSec1 == Kd::Intersection::INTERSECT_EDGE)
//      {
//         return (iSec2&1);
//      }
//      else if( iSec2 == Kd::Intersection::INTERSECT_EDGE)
//      {
//         return (iSec1&1);
//      }
//      else
//      {
//         if((iSec1&1) != (iSec2&1))
//         {
//            UBLOG(logINFO, "GbTriFaceMesh3D.isPointInGbObject3D.iSec1&1 != iSec2&1");
//            double eps = UbMath::getEqualityEpsilon<float>()*1000.0;
//            if (counter>100) {return(iSec1&1);  UBLOG(logINFO, "NACH 100 Iterationen Eps umsetzen aufgegeben!");}
//            return this->isPointInGbObject3D(x1+eps, x2+eps, x3+eps,(counter+1));
//         }
//         return iSec1&1;
//      }
//      //throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
//
//   }
//   return false;
//}
bool GbTriFaceMesh3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3, int counter)
{

    if (!nodes->empty()) {
        // Baum erstellen, wen noch keiner vorhanden
        if (!kdTree) {
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree start");
            vf::basics::Timer timer;
            timer.start();
            if (kdtreeSplitAlg == KDTREE_SAHPLIT) {
                UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit");
                this->kdTree = new Kd::Tree<double>(*this, Kd::SAHSplit<double>());
            } else if (kdtreeSplitAlg == KDTREE_SPATIALSPLIT) {
                UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SpatialMedianSplit");
                this->kdTree = new Kd::Tree<double>(*this, Kd::SpatialMedianSplit<double>());
            } else
                throw UbException(UB_EXARGS, "unknown kdtree split option)");
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - built kdTree in " << timer.getCurrentRuntimeInSeconds() << "seconds");
        }

        // eigentlicher PIO-Test
        // int iSec;
        // for(int i=0; i<100; i++)
        //{
        //   Kd::Ray<double> ray(  x1, x2, x3  //, 1, 0 ,0 );
        //                        , ( x1 < x1center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
        //                        , ( x2 < x2center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
        //                        , ( x3 < x3center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
        //                        );
        //
        //   iSec = kdTree->intersectRay( ray, Kd::CountRayIntersectionHandler<double>() );
        //
        //   if( iSec != Kd::Intersection::INTERSECT_EDGE ) //KEINE Kante getroffen
        //   {
        //      if(iSec == Kd::Intersection::ON_BOUNDARY )
        //      {
        //         return true;
        //      }
        //      return (iSec&1);  //ungerade anzahl an schnitten --> drinnen
        //   }
        //   UBLOG(logDEBUG3, "GbTriFaceMesh3D.isPointInGbObject3D.if  - an edge was hit ");
        //}
        // throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
        int iSec1, iSec2;
        double eps = 0.05;
        Kd::Ray<double> ray1(x1, x2, x3, 1.0 + eps * ((double)counter), eps * ((double)counter),
                             eps * ((double)counter));
        iSec1 = kdTree->intersectRay(ray1, Kd::CountRayIntersectionHandler<double>());
        Kd::Ray<double> ray2(x1, x2, x3, -1.0 - eps * ((double)counter), -eps * ((double)counter),
                             -eps * ((double)counter));

        iSec2 = kdTree->intersectRay(ray2, Kd::CountRayIntersectionHandler<double>());

        if (iSec1 == Kd::Intersection::ON_BOUNDARY || iSec2 == Kd::Intersection::ON_BOUNDARY) {
            return true;
        }
        if (iSec1 == Kd::Intersection::INTERSECT_EDGE && iSec2 == Kd::Intersection::INTERSECT_EDGE) {
            // UBLOG(logINFO, "GbTriFaceMesh3D.isPointInGbObject3D.INTERSECT_EDGE");

            if (counter > 20) {
                return (iSec1 & 1); /*UBLOG(logINFO, "NACH 100 Iterationen Eps umsetzen aufgegeben!");*/
            }
            return this->isPointInGbObject3D(x1, x2, x3, (counter + 1));
        } else if (iSec1 == Kd::Intersection::INTERSECT_EDGE) {
            return (iSec2 & 1);
        } else if (iSec2 == Kd::Intersection::INTERSECT_EDGE) {
            return (iSec1 & 1);
        } else {
            if ((iSec1 & 1) != (iSec2 & 1)) {
                // UBLOG(logINFO, "GbTriFaceMesh3D.isPointInGbObject3D.iSec1&1 != iSec2&1");

                if (counter > 20) {
                    return (iSec1 & 1); /* UBLOG(logINFO, "NACH 100 Iterationen Eps umsetzen aufgegeben!");*/
                }
                return this->isPointInGbObject3D(x1, x2, x3, (counter + 1));
            }
            return iSec1 & 1;
        }
        // throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
    }
    return false;
}
/*======================================================================*/
bool GbTriFaceMesh3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3)
{
    if (!nodes->empty()) {
        // Baum erstellen, wen noch keiner vorhanden
        if (!kdTree) {
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree start");
            vf::basics::Timer timer;
            timer.start();
            if (kdtreeSplitAlg == KDTREE_SAHPLIT) {
                UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit");
                //cout << "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit" << std::endl;
                this->kdTree = new Kd::Tree<double>(*this, Kd::SAHSplit<double>());
            } else if (kdtreeSplitAlg == KDTREE_SPATIALSPLIT) {
                UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SpatialMedianSplit");
                this->kdTree = new Kd::Tree<double>(*this, Kd::SpatialMedianSplit<double>());
            } else
                throw UbException(UB_EXARGS, "unknown kdtree split option)");
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - built kdTree in " << timer.getCurrentRuntimeInSeconds() << "seconds");
            //cout << "GbTriFaceMesh3D::calculateValues - built kdTree in " << timer.stop() << "seconds" << std::endl;
        }

        // eigentlicher PIO-Test
        int iSec;
        for (int i = 0; i < MAX_ITER; i++) {
            Kd::Ray<double> ray(x1, x2, x3 //, 1, 0 ,0 );
                                ,
                                (x1 < x1center ? UbRandom::rand(-1.0, -0.001, 10) : UbRandom::rand(0.001, 1.0, 10)),
                                (x2 < x2center ? UbRandom::rand(-1.0, -0.001, 10) : UbRandom::rand(0.001, 1.0, 10)),
                                (x3 < x3center ? UbRandom::rand(-1.0, -0.001, 10) : UbRandom::rand(0.001, 1.0, 10)));

            iSec = kdTree->intersectRay(ray, Kd::CountRayIntersectionHandler<double>());

            if (iSec != Kd::Intersection::INTERSECT_EDGE) // KEINE Kante getroffen
            {
                if (iSec == Kd::Intersection::ON_BOUNDARY) {
                    return true;
                }
                return (iSec & 1); // ungerade anzahl an schnitten --> drinnen
            }
            UBLOG(logDEBUG3, "GbTriFaceMesh3D.isPointInGbObject3D.if  - an edge was hit ");
        }
        throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");

        //   int iSec1,iSec2;
        //
        //   Kd::Ray<double> ray1(  x1, x2, x3, 1.0, 0.0 ,0.0 );
        //   iSec1 = kdTree->intersectRay( ray1, Kd::CountRayIntersectionHandler<double>() );
        //   Kd::Ray<double> ray2(  x1, x2, x3, -1.0, 0.0 ,0.0 );
        //   iSec2 = kdTree->intersectRay( ray2, Kd::CountRayIntersectionHandler<double>() );

        //   if(iSec1 == Kd::Intersection::ON_BOUNDARY || iSec2 == Kd::Intersection::ON_BOUNDARY)
        //   {
        //      return true;
        //   }
        //   if( iSec1 == Kd::Intersection::INTERSECT_EDGE && iSec2 == Kd::Intersection::INTERSECT_EDGE)
        //   {
        //      //UBLOG(logINFO, "GbTriFaceMesh3D.isPointInGbObject3D.INTERSECT_EDGE");
        //      double eps = UbMath::getEqualityEpsilon<double>();
        //      if (counter>100) {return(iSec1&1);  UBLOG(logINFO, "NACH 100 Iterationen Eps umsetzen aufgegeben!");}
        //      return this->isPointInGbObject3D(x1+eps, x2+eps, x3+eps,(counter+1));
        //   }
        //   else if( iSec1 == Kd::Intersection::INTERSECT_EDGE)
        //   {
        //      return (iSec2&1);
        //   }
        //   else if( iSec2 == Kd::Intersection::INTERSECT_EDGE)
        //   {
        //      return (iSec1&1);
        //   }
        //   else
        //   {
        //      if((iSec1&1) != (iSec2&1))
        //      {
        //         UBLOG(logINFO, "GbTriFaceMesh3D.isPointInGbObject3D.iSec1&1 != iSec2&1");
        //         double eps = UbMath::getEqualityEpsilon<double>();
        //         if (counter>100) {return(iSec1&1);  UBLOG(logINFO, "NACH 100 Iterationen Eps umsetzen aufgegeben!");}
        //         return this->isPointInGbObject3D(x1+eps, x2+eps, x3+eps,(counter+1));
        //      }
        //      return iSec1&1;
        //   }
        //   //throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
    }
    return false;
}
/*======================================================================*/
bool GbTriFaceMesh3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary)
{
    if (!nodes->empty()) {
        // Baum erstellen, wen noch keiner vorhanden
        if (!kdTree) {
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree start");
            vf::basics::Timer timer;
            timer.start();
            if (kdtreeSplitAlg == KDTREE_SAHPLIT) {
                UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit");
                this->kdTree = new Kd::Tree<double>(*this, Kd::SAHSplit<double>());
            } else if (kdtreeSplitAlg == KDTREE_SPATIALSPLIT) {
                UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SpatialMedianSplit");
                this->kdTree = new Kd::Tree<double>(*this, Kd::SpatialMedianSplit<double>());
            } else
                throw UbException(UB_EXARGS, "unknown kdtree split option)");
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - built kdTree in " << timer.getCurrentRuntimeInSeconds() << "seconds");
        }

        // eigentlicher PIO-Test
        int iSec;
        for (int i = 0; i < MAX_ITER; i++) {
            Kd::Ray<double> ray(
                x1, x2, x3, float((x1 < x1center ? UbRandom::rand(-1.0, -0.001, 10) : UbRandom::rand(0.001, 1.0, 10))),
                float((x2 < x2center ? UbRandom::rand(-1.0, -0.001, 10) : UbRandom::rand(0.001, 1.0, 10))),
                float((x3 < x3center ? UbRandom::rand(-1.0, -0.001, 10) : UbRandom::rand(0.001, 1.0, 10))));

            iSec = kdTree->intersectRay(ray, Kd::CountRayIntersectionHandler<double>());

            if (iSec != Kd::Intersection::INTERSECT_EDGE) // KEINE Kante getroffen
            {
                if (iSec == Kd::Intersection::ON_BOUNDARY) {
                    pointIsOnBoundary = true;
                    return true;
                }
                pointIsOnBoundary = false;
                return (iSec & 1); // ungerade anzahl an schnitten --> drinnen
            }
        }

        throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
    }

    return false;
}
/*======================================================================*/
bool GbTriFaceMesh3D::intersectLine(const double &p1_x1, const double &p1_x2, const double &p1_x3, const double &p2_x1,
                                    const double &p2_x2, const double &p2_x3)
{
    // Baum erstellen, wen noch keiner vorhanden
    if (!kdTree) {
        UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree start");
        vf::basics::Timer timer;
        timer.start();
        if (kdtreeSplitAlg == KDTREE_SAHPLIT) {
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit");
            this->kdTree = new Kd::Tree<double>(*this, Kd::SAHSplit<double>());
        } else if (kdtreeSplitAlg == KDTREE_SPATIALSPLIT) {
            UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - build KdTree with SpatialMedianSplit");
            this->kdTree = new Kd::Tree<double>(*this, Kd::SpatialMedianSplit<double>());
        } else
            throw UbException(UB_EXARGS, "unknown kdtree split option)");
        UBLOG(logDEBUG3, "GbTriFaceMesh3D::calculateValues - built kdTree in " << timer.getCurrentRuntimeInSeconds() << "seconds");
    }

    int iSec = kdTree->intersectLine(UbTupleDouble3(p1_x1, p1_x2, p1_x3), UbTupleDouble3(p2_x1, p2_x2, p2_x3),
                                     Kd::CountLineIntersectionHandler<double>());

    return (iSec != Kd::Intersection::NO_INTERSECTION);
}
/*======================================================================*/
GbLine3D *GbTriFaceMesh3D::createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
}

/*======================================================================*/
UbTuple<string, string> GbTriFaceMesh3D::writeMesh(string filename, WbWriter *writer, bool writeNormals,
                                                   vector<string> *datanames,
                                                   std::vector<std::vector<double>> *nodedata)
{
    UBLOG(logINFO, "GbTriFaceMesh3D::writeMesh ");

    vector<UbTupleFloat3> triNodes(nodes->size());
    vector<UbTupleInt3> tris(triangles->size());

    for (size_t i = 0; i < nodes->size(); i++)
        triNodes[i] = makeUbTuple((*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z);

    for (size_t i = 0; i < triangles->size(); i++)
        tris[i] = makeUbTuple((*triangles)[i].v1, (*triangles)[i].v2, (*triangles)[i].v3);

    UbTuple<string, string> filenames("", "");

    if (!datanames || datanames->empty() || !nodedata) {
        val<1>(filenames) = writer->writeTriangles(filename, triNodes, tris);
    } else {
        val<1>(filenames) = writer->writeTrianglesWithNodeData(filename, triNodes, tris, *datanames, *nodedata);
    }

    if (writeNormals) {
        vector<UbTupleFloat3> lineNodes(triangles->size() * 2);
        vector<UbTupleInt2> lines(triangles->size());
        for (size_t i = 0; i < triangles->size(); i++) {
            TriFace &triangle = (*triangles)[i];
            lineNodes[i * 2]  = makeUbTuple(triangle.getX1Centroid(*nodes), triangle.getX2Centroid(*nodes),
                                           triangle.getX3Centroid(*nodes));

            lineNodes[i * 2 + 1] = makeUbTuple((float)(triangle.getX1Centroid(*nodes) + 1.0 * triangle.nx),
                                               (float)(triangle.getX2Centroid(*nodes) + 1.0 * triangle.ny),
                                               (float)(triangle.getX3Centroid(*nodes) + 1.0 * triangle.nz));

            lines[i] = makeUbTuple((int)i * 2, (int)i * 2 + 1);
        }
        val<2>(filenames) = writer->writeLines(filename + "_normals", lineNodes, lines);
    }

    return filenames;
}
/*======================================================================*/
void GbTriFaceMesh3D::writeMeshPly(const std::string &filename)
{
    ofstream out(filename.c_str());
    if (!out)
        throw UbException(UB_EXARGS, "couldn't open " + filename);

    out << "ply" << endl;
    out << "format ascii 1.0" << endl;
    out << "element vertex " << (int)nodes->size() << endl;
    out << "property float x" << endl;
    out << "property float y" << endl;
    out << "property float z" << endl;
    out << "element face " << (int)triangles->size() << endl;
    out << "property list uchar int vertex_indices" << endl;
    out << "end_header" << endl;

    for (size_t i = 0; i < nodes->size(); i++)
        out << (*nodes)[i].x << " " << (*nodes)[i].y << " " << (*nodes)[i].z << endl;

    for (size_t i = 0; i < triangles->size(); i++)
        out << "3 " << (*triangles)[i].v1 << " " << (*triangles)[i].v2 << " " << (*triangles)[i].v3 << endl;
}
/*======================================================================*/
void GbTriFaceMesh3D::readMeshFromSTLFileASCII(string filename, bool removeRedundantNodes)
{
    UBLOG(logDEBUG1, "GbTriFaceMesh3DCreator::readMeshFromSTLFile !!! Dieses Format hat leider redundante Knoten ...");

    int nr = 0;

    ifstream in(filename.c_str());
    if (!in.good()) {
        (*nodes).clear();
        (*triangles).clear();
        UB_THROW(UbException(UB_EXARGS, "Can not open STL file: " + filename));
    }
    char title[80];
    std::string s0, s1;
    float n0, n1, n2, f0, f1, f2, f3, f4, f5, f6, f7, f8;
    in.read(title, 80);
    while (!in.eof()) {
        in >> s0; // facet || endsolid
        if (s0 == "facet") {
            in >> s1 >> n0 >> n1 >> n2; // normal x y z
            in >> s0 >> s1;             // outer loop
            in >> s0 >> f0 >> f1 >> f2; // vertex x y z
            in >> s0 >> f3 >> f4 >> f5; // vertex x y z
            in >> s0 >> f6 >> f7 >> f8; // vertex x y z
            in >> s0;                   // endloop
            in >> s0;                   // endfacet
            // Generate a new Triangle without Normal as 3 Vertices
            nodes->push_back(GbTriFaceMesh3D::Vertex(f0, f1, f2));
            nodes->push_back(GbTriFaceMesh3D::Vertex(f3, f4, f5));
            nodes->push_back(GbTriFaceMesh3D::Vertex(f6, f7, f8));
            triangles->push_back(GbTriFaceMesh3D::TriFace(nr, nr + 1, nr + 2));
            nr += 3;
        } else if (s0 == "endsolid") {
            break;
        }
    }
    in.close();

    if (removeRedundantNodes) {
        this->deleteRedundantNodes(); // dort wird autoamtisch calculateValues() aufgerufen
    } else {
        this->calculateValues();
    }
    //UBLOG(logDEBUG1, "GbTriFaceMesh3DCreator::readMeshFromSTLFile !!! Dieses Format hat leider redundante Knoten ...");

    //string dummy;

    //double x, y, z;
    //int nr = 0;

    //UbFileInputASCII in(filename);
    //in.readLine();
    //while (dummy != "endsolid") {
    //    in.readLine();
    //    in.readLine();
    //    dummy = in.readString();
    //    if (dummy != "vertex")
    //        throw UbException(UB_EXARGS, "no vertex format");
    //    x = in.readDouble();
    //    y = in.readDouble();
    //    z = in.readDouble();
    //    nodes->push_back(GbTriFaceMesh3D::Vertex((float)x, (float)y, (float)z));
    //    in.readLine();
    //    in.readString();
    //    x = in.readDouble();
    //    y = in.readDouble();
    //    z = in.readDouble();
    //    nodes->push_back(GbTriFaceMesh3D::Vertex((float)x, (float)y, (float)z));
    //    in.readLine();
    //    in.readString();
    //    x = in.readDouble();
    //    y = in.readDouble();
    //    z = in.readDouble();
    //    nodes->push_back(GbTriFaceMesh3D::Vertex((float)x, (float)y, (float)z));
    //    triangles->push_back(GbTriFaceMesh3D::TriFace(nr, nr + 1, nr + 2));
    //    in.readLine();
    //    in.readLine();
    //    in.readLine();
    //    dummy = in.readString();
    //    nr += 3;
    //    // std::cout<<"read mesh "<< nr <<" \n";
    //}

    //if (removeRedundantNodes) {
    //    this->deleteRedundantNodes(); // dort wird autoamtisch calculateValues() aufgerufen
    //} else {
    //    this->calculateValues();
    //}
}
/*======================================================================*/
void GbTriFaceMesh3D::readMeshFromSTLFileBinary(string filename, bool removeRedundantNodes)
{
    int nr  = 0;
    FILE *f = fopen(filename.c_str(), "rb");
    if (!f) {
        (*nodes).clear();
        (*triangles).clear();
        UB_THROW(UbException(UB_EXARGS, "Can not open STL file: " + filename));
    }
    char title[80];
    int nFaces;
    size_t sizef = fread(title, 80, 1, f);
    sizef        = fread((void *)&nFaces, 4, 1, f);
    float v[12]; // normal=3, vertices=3*3 = 12
    unsigned short uint16;
    // Every Face is 50 Bytes: Normal(3*float), Vertices(9*float), 2 Bytes Spacer
    for (int i = 0; i < nFaces; ++i) {
        for (size_t j = 0; j < 12; ++j) {
            sizef = fread((void *)&v[j], sizeof(float), 1, f);
        }
        sizef = fread((void *)&uint16, sizeof(unsigned short), 1, f); // spacer between successive faces
        nodes->push_back(GbTriFaceMesh3D::Vertex(v[3], v[4], v[5]));
        nodes->push_back(GbTriFaceMesh3D::Vertex(v[6], v[7], v[8]));
        nodes->push_back(GbTriFaceMesh3D::Vertex(v[9], v[10], v[11]));
        triangles->push_back(GbTriFaceMesh3D::TriFace(nr, nr + 1, nr + 2));
        nr += 3;
    }
    (void)sizef;
    fclose(f);

    if (removeRedundantNodes) {
        this->deleteRedundantNodes(); // dort wird autoamtisch calculateValues() aufgerufen
    } else {
        this->calculateValues();
    }
}

//! \}
