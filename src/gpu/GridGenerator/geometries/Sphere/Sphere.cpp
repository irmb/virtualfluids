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
//! \addtogroup gpu_geometries geometries
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "Sphere.h"

#include <algorithm>    // std::min
#include <float.h>
#include <cmath>

#include "geometries/Vertex/Vertex.h"

Sphere::Sphere(const double& centerX, const double& centerY, const double& centerZ, const double& radius)
    : centerX(centerX), centerY(centerY), centerZ(centerZ), radius(radius)
{

}

SPtr<Sphere> Sphere::makeShared(double centerX, double centerY, double centerZ, double radius)
{
    return std::make_shared<Sphere>(centerX, centerY, centerZ, radius);
}

SPtr<Object> Sphere::clone() const
{
    return std::make_shared<Sphere>(centerX, centerY, centerZ, radius);
}

double Sphere::getX1Centroid() const
{
    return centerX;
}

double Sphere::getX1Minimum() const
{
    return centerX - radius;
}

double Sphere::getX1Maximum() const
{
    return centerX + radius;
}

double Sphere::getX2Centroid() const
{
    return centerY;
}

double Sphere::getX2Minimum() const
{
    return centerY - radius;
}

double Sphere::getX2Maximum() const
{
    return centerY + radius;
}

double Sphere::getX3Centroid() const
{
    return centerZ;
}

double Sphere::getX3Minimum() const
{
    return centerZ - radius;
}

double Sphere::getX3Maximum() const
{
    return centerZ + radius;
}

bool Sphere::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset,
    const double& maxOffset)
{
    double offset = maxOffset;
    if (x1 < centerX || x2 < centerY || x3 < centerZ)
        offset = minOffset;
        

    const double deltaX1 = x1 - centerX;
    const double deltaX2 = x2 - centerY;
    const double deltaX3 = x3 - centerZ;

    return (deltaX1*deltaX1 + deltaX2*deltaX2 + deltaX3*deltaX3) < ((this->radius - offset) * (this->radius - offset));
}


void Sphere::changeSizeByDelta(double delta)
{
    this->radius += delta;
}

int Sphere::getIntersection(const Vertex & point, const Vertex & direction, Vertex & pointOnObject, real & qVal)
{
    
    Vertex relativePoint( point.x - this->centerX, 
                          point.y - this->centerY, 
                          point.z - this->centerZ );

    real directionSquare = direction.x * direction.x
                         + direction.y * direction.y
                         + direction.z * direction.z;

    real p = 2* ( relativePoint.x * direction.x
                + relativePoint.y * direction.y
                + relativePoint.z * direction.z )
                / directionSquare;

    real q = ( relativePoint.x * relativePoint.x
             + relativePoint.y * relativePoint.y
             + relativePoint.z * relativePoint.z
             - (real)this->radius * (real)this->radius )
           / directionSquare;

    real discriminant = 0.25 * p * p - q;


    if( discriminant < 0.0 ) return 1;

    real result1 = - 0.5 * p + std::sqrt( discriminant );
    real result2 = - 0.5 * p - std::sqrt( discriminant );

    if( result1 < 0.0 && result2 < 0.0 ) return 1;

    if (result1 < 0.0)
        result1 = (real)FLT_MAX;
    if (result2 < 0.0)
        result2 = (real)FLT_MAX;

    real t = std::min( result1, result2 );

    pointOnObject.x = point.x + t * direction.x;
    pointOnObject.y = point.y + t * direction.y;
    pointOnObject.z = point.z + t * direction.z;

    qVal = t;

    return 0;
}

//! \}
