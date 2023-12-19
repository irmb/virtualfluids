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
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include <GbPoint3D.h>
using namespace std;

/*=======================================================*/
GbPoint3D::GbPoint3D()
{
    this->x1 = 0.0;
    this->x2 = 0.0;
    this->x3 = 0.0;
}
/*=======================================================*/
GbPoint3D::GbPoint3D(const double &x1, const double &x2, const double &x3)
{
    this->x1 = x1;
    this->x2 = x2;
    this->x3 = x3;
}
/*=======================================================*/
GbPoint3D::GbPoint3D(GbPoint3D *point)
{
    this->x1 = point->x1;
    this->x2 = point->x2;
    this->x3 = point->x3;
}
/*=======================================================*/
double GbPoint3D::getDistance(GbPoint3D *p)
{
    double dx1 = this->x1 - p->x1;
    double dx2 = this->x2 - p->x2;
    double dx3 = this->x3 - p->x3;
    return std::sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
}
/*=======================================================*/
bool GbPoint3D::equals(const GbPoint3D *point) const
{
    if (fabs(this->x1 - point->x1) > 1.E-10)
        return false;
    if (fabs(this->x2 - point->x2) > 1.E-10)
        return false;
    if (fabs(this->x3 - point->x3) > 1.E-10)
        return false;

    return true;
}
/*=======================================================*/
void GbPoint3D::transform(const double matrix[4][4])
{
    double tempX1 = x1;
    double tempX2 = x2;
    double tempX3 = x3;
    x1            = matrix[0][0] * tempX1 + matrix[0][1] * tempX2 + matrix[0][2] * tempX3 + matrix[0][3] * 1.;
    x2            = matrix[1][0] * tempX1 + matrix[1][1] * tempX2 + matrix[1][2] * tempX3 + matrix[1][3] * 1.;
    x3            = matrix[2][0] * tempX1 + matrix[2][1] * tempX2 + matrix[2][2] * tempX3 + matrix[2][3] * 1.;
}
/*=======================================================*/
bool GbPoint3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3)
{
    return (fabs(x1) < 1.E-13 && fabs(x2) < 1.E-13 && fabs(x3) < 1.E-13);
}
/*=======================================================*/
bool GbPoint3D::isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary)
{
    pointIsOnBoundary = (fabs(x1) < 1.E-13 && fabs(x2) < 1.E-13 && fabs(x3) < 1.E-13);
    return pointIsOnBoundary;
}
/*=======================================================*/
vector<GbTriangle3D *> GbPoint3D::getSurfaceTriangleSet()
{
    cout << "GbPoint3D::getSurfaceTriangleSet() - test ... if no exception occurs, everything is fine\n";
    vector<GbTriangle3D *> triangles;
    return triangles; //<-empty vector! is okay!

    // old:
    // to avoid unnecessary exceptions a point will generate a triangle with
    // c3 point with same coordinates
    // vector<GbTriangle3D*> triangles;
    // GbPoint3D p1(getX1Coordinate(),getX2Coordinate(),getX3Coordinate());
    // triangles.push_back(new GbTriangle3D(new GbPoint3D(p1),new GbPoint3D(p1),new GbPoint3D(p1)));
}
/*=======================================================*/
GbLine3D *GbPoint3D::createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/)
{
    throw UbException(UB_EXARGS, "not implemented");
}
/*=======================================================*/
string GbPoint3D::toString()
{
    stringstream ss;
    ss << "GbPoint3D[" << this->x1 << "," << this->x2 << "," << this->x3 << "]";
    return ((ss.str()).c_str());
}
/*=======================================================*/
void GbPoint3D::translate(const double &dx1, const double &dx2, const double &dx3)
{
    this->x1 += dx1;
    this->x2 += dx2;
    this->x3 += dx3;
    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbPoint3D::rotate(const double &rx1, const double &rx2, const double &rx3)
{
    double newX1 = cos(rx3) * cos(rx2) * x1 - x2 * sin(rx3) * cos(rx1) + x2 * cos(rx3) * sin(rx2) * sin(rx1) +
                   x3 * sin(rx3) * sin(rx1) + x3 * cos(rx3) * sin(rx2) * cos(rx1);
    double newX2 = sin(rx3) * cos(rx2) * x1 + x2 * cos(rx3) * cos(rx1) + x2 * sin(rx3) * sin(rx2) * sin(rx1) -
                   x3 * cos(rx3) * sin(rx1) + x3 * sin(rx3) * sin(rx2) * cos(rx1);
    double newX3 = -sin(rx2) * x1 + cos(rx2) * sin(rx1) * x2 + cos(rx2) * cos(rx1) * x3;

    this->x1 = newX1;
    this->x2 = newX2;
    this->x3 = newX3;
    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbPoint3D::scale(const double &sx1, const double &sx2, const double &sx3)
{
    this->x1 *= sx1;
    this->x2 *= sx2;
    this->x3 *= sx3;
    this->notifyObserversObjectChanged();
}
/*=======================================================*/

//! \}
