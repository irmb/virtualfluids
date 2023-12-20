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
#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbObjectGroup3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>

using namespace std;

/*=====================================================*/
GbObjectGroup3D::GbObjectGroup3D() { this->setName("ObjectGroup"); }
/*=====================================================*/
GbObjectGroup3D::~GbObjectGroup3D() = default;
/*=====================================================*/
void GbObjectGroup3D::finalize() { throw UbException(UB_EXARGS, "not implemented."); }
/*=======================================================*/
void GbObjectGroup3D::setCenterCoordinates(const double & /*x1*/, const double & /*x2*/, const double & /*x3*/)
{
    throw UbException(UB_EXARGS, "not implemented.");
}
/*=====================================================*/
double GbObjectGroup3D::getDistance(GbPoint3D * /*p*/) { throw UbException(UB_EXARGS, "not implemented."); }
/*=====================================================*/

void GbObjectGroup3D::setCenterX1Coordinate(const double & /*value*/)
{
    throw UbException(UB_EXARGS, "not implemented.");
}
/*=====================================================*/
void GbObjectGroup3D::setCenterX2Coordinate(const double & /*value*/)
{
    throw UbException(UB_EXARGS, "not implemented.");
}
/*=====================================================*/
void GbObjectGroup3D::setCenterX3Coordinate(const double & /*value*/)
{
    throw UbException(UB_EXARGS, "not implemented.");
}
/*=====================================================*/
void GbObjectGroup3D::setRadius(const double & /*radius*/) { throw UbException(UB_EXARGS, "not implemented."); }
/*=====================================================*/
double GbObjectGroup3D::getDistance(const double & /*x1p*/, const double & /*x2p*/, const double & /*x3p*/)
{
    throw UbException(UB_EXARGS, "not implemented.");
}
/*=====================================================*/
// true, wenn 'in Object' oder 'auf Boundary'!
bool GbObjectGroup3D::isPointInGbObject3D(const double & /*x1p*/, const double & /*x2p*/, const double & /*x3p*/)
{
    return false;
}
/*=====================================================*/
// true, wenn 'in Object' oder 'auf Boundary'!
bool GbObjectGroup3D::isPointInGbObject3D(const double & /*x1p*/, const double & /*x2p*/, const double & /*x3p*/,
                                          bool & /*pointIsOnBoundary*/)
{
    return false;
}
/*=====================================================*/
string GbObjectGroup3D::toString()
{
    stringstream ss;
    ss << "GbObjectGroup3D[";
    ss << "mid=" << midPoint->toString() << ", r=" << radius << "]";
    return ss.str();
}
/*=====================================================*/
GbLine3D *GbObjectGroup3D::createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/) { return NULL; }
/*=========================================================================*/
vector<GbTriangle3D *> GbObjectGroup3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> allTriangles;

    // loop ueber alle objekte in der group
    for (std::list<GbObject3D *>::iterator iter = this->geoobjects.begin(); iter != this->geoobjects.end(); iter++) {
        vector<GbTriangle3D *> triangles;
        triangles = (*iter)->getSurfaceTriangleSet();

        for (size_t i = 0; i < triangles.size(); i++) {
            // kopieren...
            allTriangles.push_back(triangles[i]);
        }
    }
    return allTriangles;
}
/*=======================================================*/
void GbObjectGroup3D::addSurfaceTriangleSet(vector<UbTupleFloat3> &/*nodes*/, vector<UbTupleInt3> &/*triangles*/) {}
/*=======================================================*/
bool GbObjectGroup3D::hasIntersectionWithDirectedLine(GbPoint3D /*origin*/, GbPoint3D /*direction*/) { return false; }
/*=======================================================*/
bool GbObjectGroup3D::isCellCuttingGbObject3D(const double & /*x1a*/, const double & /*x2a*/, const double & /*x3a*/,
                                              const double & /*x1b*/, const double & /*x2b*/, const double & /*x3b*/)
{
    return false;
}
/*=======================================================*/
bool GbObjectGroup3D::isCellInsideOrCuttingGbObject3D(const double & /*x1a*/, const double & /*x2a*/,
                                                      const double & /*x3a*/, const double & /*x1b*/,
                                                      const double & /*x2b*/, const double & /*x3b*/)
{
    return false;
}
/*==========================================================*/
double GbObjectGroup3D::getCellVolumeInsideGbObject3D(const double & /*x1a*/, const double & /*x2a*/,
                                                      const double & /*x3a*/, const double & /*x1b*/,
                                                      const double & /*x2b*/, const double & /*x3b*/)
{
    return 0.0;
}
/*==========================================================*/
double GbObjectGroup3D::getIntersectionRaytraceFactor(const double & /*x1*/, const double & /*x2*/,
                                                      const double & /*x3*/, const double & /*rx1*/,
                                                      const double & /*rx2*/, const double & /*rx3*/)
{
    return 0.0;
}
/*=======================================================*/

//! \}
