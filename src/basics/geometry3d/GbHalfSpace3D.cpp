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
#include <geometry3d/GbHalfSpace3D.h>

using namespace std;

/*==========================================================*/
GbHalfSpace3D::GbHalfSpace3D(GbTriangle3D *triangle)
{
    GbPoint3D *PointA = triangle->getPoint1();
    GbPoint3D *PointB = triangle->getPoint2();
    GbPoint3D *PointC = triangle->getPoint3();

    GbVector3D A(PointA->x1, PointA->x2, PointA->x3);
    GbVector3D BA(PointB->x1 - PointA->x1, PointB->x2 - PointA->x2, PointB->x3 - PointA->x3);
    GbVector3D CA(PointC->x1 - PointA->x1, PointC->x2 - PointA->x2, PointC->x3 - PointA->x3);
    GbVector3D BACA = BA.Cross(CA);
    // this->Normal = PointB->subtract(PointA)->cross(PointC->subtract(PointA))->normalize();
    BACA.Normalize();
    // this->Normal = BACA;
    normalX = BACA[0];
    normalY = BACA[1];
    normalZ = BACA[2];
    // this->d = this->Normal.Dot(A);
    this->d = normalX * A[0] + normalY * A[1] + normalZ * A[2];
}
/*==========================================================*/
GbHalfSpace3D::GbHalfSpace3D(GbPoint3D *PointA, GbPoint3D *PointB, GbPoint3D *PointC)
{
    GbVector3D A(PointA->x1, PointA->x2, PointA->x3);
    GbVector3D BA(PointB->x1 - PointA->x1, PointB->x2 - PointA->x2, PointB->x3 - PointA->x3);
    GbVector3D CA(PointC->x1 - PointA->x1, PointC->x2 - PointA->x2, PointC->x3 - PointA->x3);
    GbVector3D BACA = BA.Cross(CA);
    // this->Normal = PointB->subtract(PointA)->cross(PointC->subtract(PointA))->normalize();
    BACA.Normalize();
    // this->Normal = BACA;
    normalX = BACA[0];
    normalY = BACA[1];
    normalZ = BACA[2];
    // this->d = this->Normal.Dot(A);
    this->d = normalX * A[0] + normalY * A[1] + normalZ * A[2];
}
/*==========================================================*/
GbHalfSpace3D::GbHalfSpace3D(GbPoint3D *PointA, GbPoint3D *PointB)
{
    GbVector3D A(PointA->x1, PointA->x2, PointA->x3);
    GbVector3D B(PointB->x1, PointB->x2, PointB->x3);
    GbVector3D K(0.0, 0.0, 0.99); // the vector from PointA - third point

    GbVector3D PointBA  = B - A;
    GbVector3D PointBAK = PointBA.Cross(K);
    PointBAK.Normalize();

    // this->Normal = PointBAK;
    normalX = PointBAK[0];
    normalY = PointBAK[1];
    normalZ = PointBAK[2];

    // this->d = this->Normal.Dot(A);
    this->d = normalX * PointA->x1 + normalY * PointA->x2 + normalZ * PointA->x3;
}
/*==========================================================*/
GbHalfSpace3D::GbHalfSpace3D(const double &p1x, const double &p1y, const double &p1z, const double &p2x,
                             const double &p2y, const double &p2z, const double &p3x, const double &p3y,
                             const double &p3z)
{
    double p2minusP1x = p2x - p1x;
    double p2minusP1y = p2y - p1y;
    double p2minusP1z = p2z - p1z;

    double P3minusP1x = p3x - p1x;
    double P3minusP1y = p3y - p1y;
    double P3minusP1z = p3z - p1z;

    // normal = BA x CA
    normalX = p2minusP1y * P3minusP1z - p2minusP1z * P3minusP1y;
    normalY = p2minusP1z * P3minusP1x - p2minusP1x * P3minusP1z;
    normalZ = p2minusP1x * P3minusP1y - p2minusP1y * P3minusP1x;

    // normalize BACA
    double oneOverNormalLength = 1.0 / (std::sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ));
    normalX *= oneOverNormalLength;
    normalY *= oneOverNormalLength;
    normalZ *= oneOverNormalLength;

    // d = normal * p1
    this->d = normalX * p1x + normalY * p1y + normalZ * p1z;
}
/*==========================================================*/
GbHalfSpace3D::GbHalfSpace3D(const double &p1x, const double &p1y, const double &p1z, const double &nx,
                             const double &ny, const double &nz)
{
    // normal = BA x CA
    normalX = nx;
    normalY = ny;
    normalZ = nz;

    // d = normal * p1
    this->d = nx * p1x + ny * p1y + nz * p1z;
}
/*==========================================================*/

//! \}
