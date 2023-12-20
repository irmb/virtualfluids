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
#include "VerticalCylinder.h"

VerticalCylinder::VerticalCylinder(const double& centerX, const double& centerY, const double& centerZ, const double& radius, const double& height)
    : centerX(centerX), centerY(centerY), centerZ(centerZ), radius(radius), height(height)
{

}

SPtr<VerticalCylinder> VerticalCylinder::makeShared(double centerX, double centerY, double centerZ, double radius, double height)
{
    return std::make_shared<VerticalCylinder>(centerX, centerY, centerZ, radius, height);
}

SPtr<Object> VerticalCylinder::clone() const
{
    return std::make_shared<VerticalCylinder>(centerX, centerY, centerZ, radius, height);
}

double VerticalCylinder::getX1Centroid() const
{
    return centerX;
}

double VerticalCylinder::getX1Minimum() const
{
    return centerX - radius;
}

double VerticalCylinder::getX1Maximum() const
{
    return centerX + radius;
}

double VerticalCylinder::getX2Centroid() const
{
    return centerY;
}

double VerticalCylinder::getX2Minimum() const
{
    return centerY - radius;
}

double VerticalCylinder::getX2Maximum() const
{
    return centerY + radius;
}

double VerticalCylinder::getX3Centroid() const
{
    return centerZ;
}

double VerticalCylinder::getX3Minimum() const
{
    return centerZ - 0.5 * height;
}

double VerticalCylinder::getX3Maximum() const
{
    return centerZ + 0.5 * height;
}

bool VerticalCylinder::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
    double offset = maxOffset;
    if (x1 < centerX || x2 < centerY || x3 < centerZ)
        offset = minOffset;
        

    const double deltaX1 = x1 - centerX;
    const double deltaX2 = x2 - centerY;
    const double deltaX3 = x3 - centerZ;

    if( deltaX3 > 0.5 * height || deltaX3 < - 0.5 * height )
        return false;

    return (deltaX1*deltaX1 + deltaX2*deltaX2) < ((this->radius - offset) * (this->radius - offset));
}


void VerticalCylinder::changeSizeByDelta(double delta)
{
    this->radius += delta;
    this->height += 2 * delta;
}

double VerticalCylinder::getRadius() const
{
    return radius;
}

double VerticalCylinder::getHeight() const
{
    return height;
}
//! \}
