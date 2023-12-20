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
#include <basics/utilities/UbMath.h>
#include <geometry3d/CoordinateTransformation3D.h>

using namespace std;

CoordinateTransformation3D::CoordinateTransformation3D()
{
    this->setTransformationValues(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
}
/*======================================================*/
CoordinateTransformation3D::CoordinateTransformation3D(const double &originX1, const double &originX2,
                                                       const double &originX3, const double &dx1, const double &dx2,
                                                       const double &dx3, const double &alpha, const double &beta,
                                                       const double &gamma)
{
    this->setTransformationValues(originX1, originX2, originX3, dx1, dx2, dx3, alpha, beta, gamma);
}
/*======================================================*/
CoordinateTransformation3D::CoordinateTransformation3D(const double &originX1, const double &originX2,
                                                       const double &originX3, const double &dx1, const double &dx2,
                                                       const double &dx3)
{
    this->setTransformationValues(originX1, originX2, originX3, dx1, dx2, dx3, 0.0, 0.0, 0.0);
}
/*======================================================*/
CoordinateTransformation3D::CoordinateTransformation3D(CoordinateTransformation3D *transformation)
{
    this->setTransformationValues(transformation->Tx1, transformation->Tx2, transformation->Tx3, transformation->Sx1,
                                  transformation->Sx2, transformation->Sx3, transformation->alpha, transformation->beta,
                                  transformation->gamma);
}
/*======================================================*/
// void CoordinateTransformation3D::init()
// {
//    this->Tx1   = 0.0;      this->Tx2   = 0.0;    this->Tx3   = 0.0;
//    this->Sx1   = 1.0;      this->Sx2   = 1.0;    this->Sx3   = 1.0;
//    this->alpha = 0.0;        this->beta = 0.0;        this->gamma = 0.0;
//
//    this->toX1factorX1   = 1.0; this->toX1factorX2   = 0.0; this->toX1factorX3   = 0.0;
//    this->toX2factorX1   = 0.0; this->toX2factorX2   = 1.0; this->toX2factorX3   = 0.0;
//    this->toX3factorX1   = 0.0; this->toX3factorX2   = 0.0; this->toX3factorX3   = 1.0;
//    this->toX1delta      = 0.0; this->toX2delta      = 0.0; this->toX3delta      = 0.0;
//    this->fromX1factorX1 = 1.0; this->fromX1factorX2 = 0.0; this->fromX1factorX3 = 0.0;
//    this->fromX2factorX1 = 0.0; this->fromX2factorX2 = 1.0; this->fromX2factorX3 = 0.0;
//    this->fromX3factorX1 = 0.0; this->fromX3factorX2 = 0.0; this->fromX3factorX3 = 1.0;
//
//    this->active         = false;
//    this->transformation = false;
// }
/*======================================================*/

/**====  Set transformation values  ====**/
/*!
\brief Set transformation values
@param a     transformed coordinate system x0 (in global coordinates)
@param b     transformed coordinate system y0 (in global coordinates)
@param c     transformed coordinate system z0 (in global coordinates)
@param dx1    x coordinate scaling       (dx_transformed/dx_global)
@param dx2    y coordinate scaling       (dy_transformed/dy_global)
@param dx3    z coordinate scaling       (dz_transformed/dz_global)
@param alpha rotation around z angle    (positive FROM global TO transformed coordinate system)
@param beta  rotation around y angle
@param gamma rotation around x angle
@exception IllegalArgumentException if c1 of the scale values is between -1.0E-8 and 1.0E-8
*/

void CoordinateTransformation3D::setTransformationValues(const double &originX1, const double &originX2,
                                                         const double &originX3, const double &dx1, const double &dx2,
                                                         const double &dx3, const double &alpha, const double &beta,
                                                         const double &gamma)
{
    if (UbMath::zero(dx1) || UbMath::zero(dx2) || UbMath::zero(dx3))
        throw UbException(UB_EXARGS, "error: at least one delta==0.0");

    this->Tx1   = originX1;
    this->Tx2   = originX2;
    this->Tx3   = originX3;
    this->Sx1   = dx1;
    this->Sx2   = dx2;
    this->Sx3   = dx3;
    this->alpha = alpha;
    this->beta  = beta;
    this->gamma = gamma;

    double ra   = UbMath::PI * alpha / 180.0;
    double cosA = cos(ra);
    double sinA = sin(ra);
    double rb   = UbMath::PI * beta / 180.0;
    double cosB = cos(rb);
    double sinB = sin(rb);
    double rg   = UbMath::PI * gamma / 180.0;
    double cosG = cos(rg);
    double sinG = sin(rg);

    // Matrix-Werte von T_invers  (indizes: 12 = spalte 1 zeile 2)
    double divisor = (Sx1 * Sx2 * Sx3);

    this->toX1factorX1 = +cosB * cosA * Sx2 * Sx3 / divisor;
    this->toX1factorX2 = -cosB * sinA * Sx1 * Sx3 / divisor;
    this->toX1factorX3 = +sinB * Sx1 * Sx2 / divisor;
    this->toX1delta =
        (-Tx3 * Sx1 * Sx2 * sinB + Tx2 * Sx1 * Sx3 * sinA * cosB - Tx1 * Sx2 * Sx3 * cosB * cosA) / divisor;

    this->toX2factorX1 = Sx2 * Sx3 * (sinG * sinB * cosA + cosG * sinA) / divisor;
    this->toX2factorX2 = Sx1 * Sx3 * (-sinG * sinB * sinA + cosG * cosA) / divisor;
    this->toX2factorX3 = -Sx1 * Sx2 * cosB * sinG / divisor;
    this->toX2delta =
        (-Tx2 * Sx1 * Sx3 * cosG * cosA + Tx3 * Sx1 * Sx2 * sinG * cosB + Tx2 * Sx1 * Sx3 * sinG * sinA * sinB -
         Tx1 * Sx2 * Sx3 * cosG * sinA - Tx1 * Sx2 * Sx3 * sinB * sinG * cosA) /
        divisor;

    this->toX3factorX1 = Sx2 * Sx3 * (-cosG * sinB * cosA + sinG * sinA) / divisor;
    this->toX3factorX2 = Sx1 * Sx3 * (sinB * cosG * sinA + sinG * cosA) / divisor;
    this->toX3factorX3 = Sx1 * Sx2 * cosB * cosG / divisor;
    this->toX3delta =
        (-Tx2 * Sx1 * Sx3 * sinG * cosA - Tx3 * Sx1 * Sx2 * cosG * cosB - Tx2 * Sx1 * Sx3 * cosG * sinA * sinB -
         Tx1 * Sx2 * Sx3 * sinG * sinA + Tx1 * Sx2 * Sx3 * sinB * cosG * cosA) /
        divisor;

    // Matrix-Werte von T_invers  (indizes: 12 = spalte 1 zeile 2)
    this->fromX1factorX1 = cosB * cosA * Sx1;
    this->fromX1factorX2 = (sinG * sinB * cosA + cosG * sinA) * Sx1;
    this->fromX1factorX3 = (-cosG * sinB * cosA + sinG * sinA) * Sx1;
    this->fromX1delta    = Tx1;

    this->fromX2factorX1 = -cosB * sinA * Sx2;
    this->fromX2factorX2 = -(sinG * sinB * sinA - cosG * cosA) * Sx2;
    this->fromX2factorX3 = (cosG * sinB * sinA + sinG * cosA) * Sx2;
    this->fromX2delta    = Tx2;

    this->fromX3factorX1 = sinB * Sx3;
    this->fromX3factorX2 = -sinG * cosB * Sx3;
    this->fromX3factorX3 = cosG * cosB * Sx3;
    this->fromX3delta    = Tx3;

    this->active = true;

    this->transformation = true;
}
/*======================================================*/
/*!
Set transformation active state (if this IS a transformation)
@param active true to be active, false otherwise
**/
void CoordinateTransformation3D::setActive(const bool &active)
{
    if (this->active == active)
        return;
    if (this->transformation)
        this->active = active;
}
/*======================================================*/
/*!
Transform FROM global coordinates TO transformed coordinates.
@param x1  the global x coordinate
@param x2  the global y coordinate
@param x3  the global z coordinate
**/
double CoordinateTransformation3D::transformForwardToX1Coordinate(const double &x1, const double &x2,
                                                                  const double &x3) const
{
    if (this->active)
        return this->toX1factorX1 * x1 + this->toX1factorX2 * x2 + this->toX1factorX3 * x3 + this->toX1delta;
    else
        return x1;
}
/*======================================================*/
double CoordinateTransformation3D::transformForwardToX2Coordinate(const double &x1, const double &x2,
                                                                  const double &x3) const
{
    if (this->active)
        return this->toX2factorX1 * x1 + this->toX2factorX2 * x2 + this->toX2factorX3 * x3 + this->toX2delta;
    else
        return x2;
}
/*======================================================*/
double CoordinateTransformation3D::transformForwardToX3Coordinate(const double &x1, const double &x2,
                                                                  const double &x3) const
{
    if (this->active)
        return this->toX3factorX1 * x1 + this->toX3factorX2 * x2 + this->toX3factorX3 * x3 + this->toX3delta;
    else
        return x3;
}
/*======================================================*/
/*!
Transform FROM global coordinates TO transformed coordinates (ignoring rotation).
@param x1  the global x coordinate
**/
double CoordinateTransformation3D::transformForwardToX1CoordinateIgnoringRotation(const double &x1) const
{
    if (this->active)
        return (x1 - this->Tx1) / this->Sx1;
    else
        return x1;
}
/*======================================================*/
double CoordinateTransformation3D::transformForwardToX2CoordinateIgnoringRotation(const double &x2) const
{
    if (this->active)
        return (x2 - this->Tx2) / this->Sx2;
    else
        return x2;
}
/*======================================================*/
double CoordinateTransformation3D::transformForwardToX3CoordinateIgnoringRotation(const double &x3) const
{
    if (this->active)
        return (x3 - this->Tx3) / this->Sx3;
    else
        return x3;
}
/*======================================================*/
/*!
Transform FROM transformed coordinates TO global coordinates.
@param x1  the transformed x coordinate
@param x2  the transformed y coordinate
@param x3  the transformed z coordinate
**/
double CoordinateTransformation3D::transformBackwardToX1Coordinate(const double &x1, const double &x2,
                                                                   const double &x3) const
{
    if (this->active)
        return this->fromX1factorX1 * x1 + this->fromX1factorX2 * x2 + this->fromX1factorX3 * x3 + this->fromX1delta;
    else
        return x1;
}
/*======================================================*/
double CoordinateTransformation3D::transformBackwardToX2Coordinate(const double &x1, const double &x2,
                                                                   const double &x3) const
{
    if (this->active)
        return this->fromX2factorX1 * x1 + this->fromX2factorX2 * x2 + this->fromX2factorX3 * x3 + this->fromX2delta;
    else
        return x2;
}
/*======================================================*/
double CoordinateTransformation3D::transformBackwardToX3Coordinate(const double &x1, const double &x2,
                                                                   const double &x3) const
{
    if (this->active)
        return this->fromX3factorX1 * x1 + this->fromX3factorX2 * x2 + this->fromX3factorX3 * x3 + this->fromX3delta;
    else
        return x3;
}
/*======================================================*/
/*!
Transform FROM transformed coordinates TO global coordinates (ignoring rotation).
@param x1  the transformed x coordinate
**/
double CoordinateTransformation3D::transformBackwardToX1CoordinateIgnoringRotation(const double &x1) const
{
    if (this->active)
        return x1 * this->Sx1 + this->Tx1;
    else
        return x1;
}
/*======================================================*/
double CoordinateTransformation3D::transformBackwardToX2CoordinateIgnoringRotation(const double &x2) const
{
    if (this->active)
        return x2 * this->Sx2 + this->Tx2;
    else
        return x2;
}
/*======================================================*/
double CoordinateTransformation3D::transformBackwardToX3CoordinateIgnoringRotation(const double &x3) const
{
    if (this->active)
        return x3 * this->Sx3 + this->Tx3;
    else
        return x3;
}
/*======================================================*/
/*!
Returns a string representation of this transformation.
@return a string representation of this transformation
**/
string CoordinateTransformation3D::toString() const
{
    stringstream ss;
    ss << " CoordinateTransformation3D\n";
    //    ss<<"[isTransformation="<<this->transformation;
    //    ss<<", isActive="<<this->active<<endl;
    ss << " ,a=" << this->Tx1 << ", b=" << this->Tx2 << ", c=" << this->Tx3 << endl;
    ss << " , dx1=" << this->Sx1 << ", dx2=" << this->Sx2 << ", dx2=" << this->Sx3 << endl;
    //    ss<<" , alpha="<<this->alpha<<", beta="<<this->beta<endl;
    //    ss<<"]";
    //    ss<<"[to11="<<this->to11<<", to12="<<this->to12<<", to13="<<this->to13;
    //    ss<<", to21="<<this->to21<<", to22="<<this->to22<<", to23="<<this->to23;
    //    ss<<", to31="<<this->to31<<", to32="<<this->to32<<", to33="<<this->to33;
    //    ss<<", toA="<<this->toA<<", toB="<<this->toB<<", toC="<<this->toC;
    //    ss<<", from11="<<this->from11<<", from12="<<this->from12<<", from13="<<this->from13;
    //    ss<<", from21="<<this->from21<<", from22="<<this->from22<<", from23="<<this->from23;
    //    ss<<", from31="<<this->from31<<", from32="<<this->from32<<", from33="<<this->from33;
    //    ss<<", fromA="<<this->fromA; ss<<", fromB="<<this->fromB; ss<<", fromC="<<this->fromC;
    //    ss<<"]}";
    return ss.str();
}

//! \}
