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
//=======================================================================================
#include "Cylinder.h"
#include <numeric>

using namespace axis;

Cylinder::Cylinder(double centerX, double centerY, double centerZ, double radius, double height, Axis rotationalAxis)
    : center({ centerX, centerY, centerZ }), radius(radius), height(height), rotationalAxis(rotationalAxis)
{
}

Cylinder::Cylinder(std::array<double, 3> center, double radius, double height, Axis axis)
    : center(center), radius(radius), height(height), rotationalAxis(axis)
{
}

SPtr<Object> Cylinder::clone() const
{
    return std::make_shared<Cylinder>(center, radius, height, rotationalAxis);
}

double Cylinder::getCentroidCoordinate(Axis coordinateDirection) const
{
    return center.at(coordinateDirection);
}

double Cylinder::getMinimunCoordinate(Axis coordinateDirection) const
{
    const auto unitVector = unitVectors.at(rotationalAxis);
    return center.at(coordinateDirection) - 0.5 * height * unitVector.at(coordinateDirection) +
           radius * (unitVector.at(coordinateDirection) - 1);
}

double Cylinder::getMaximumCoordinate(Axis coordinateDirection) const
{
    const auto unitVector = unitVectors.at(rotationalAxis);
    return center.at(coordinateDirection) + 0.5 * height * unitVector.at(coordinateDirection) -
           radius * (unitVector.at(coordinateDirection) - 1);
}

double Cylinder::getX1Centroid() const
{
    return getCentroidCoordinate(x);
}

double Cylinder::getX1Minimum() const
{
    return getMinimunCoordinate(x);
}

double Cylinder::getX1Maximum() const
{
    return getMaximumCoordinate(x);
}

double Cylinder::getX2Centroid() const
{
    return getCentroidCoordinate(y);
}

double Cylinder::getX2Minimum() const
{
    return getMinimunCoordinate(y);
}

double Cylinder::getX2Maximum() const
{
    return getMaximumCoordinate(y);
}

double Cylinder::getX3Centroid() const
{
    return getCentroidCoordinate(z);
}

double Cylinder::getX3Minimum() const
{
    return getMinimunCoordinate(z);
}

double Cylinder::getX3Maximum() const
{
    return getMaximumCoordinate(z);
}

double Cylinder::getRadius() const
{
    return radius;
}

double Cylinder::getHeight() const
{
    return height;
}

Axis Cylinder::getRotationalAxis() const
{
    return rotationalAxis;
}

bool Cylinder::isInCircle(double delta1, double delta2, double offset) const
{
    return (delta1 * delta1 + delta2 * delta2) < ((this->radius - offset) * (this->radius - offset));
}

bool Cylinder::isPointInObject(const double &x1, const double &x2, const double &x3, const double &minOffset,
                               const double &maxOffset)
{
    double offset = maxOffset;
    if (x1 < center.at(x) || x2 < center.at(y) || x3 < center.at(z)) offset = minOffset;

    const double deltaX1 = x1 - center.at(x);
    const double deltaX2 = x2 - center.at(y);
    const double deltaX3 = x3 - center.at(z);

    switch (rotationalAxis) {
        case x:
            if (deltaX1 > 0.5 * height || deltaX1 < -0.5 * height) return false;
            return isInCircle(deltaX2, deltaX3, offset);
        case y:
            if (deltaX2 > 0.5 * height || deltaX2 < -0.5 * height) return false;
            return isInCircle(deltaX1, deltaX3, offset);
        case z:
            if (deltaX3 > 0.5 * height || deltaX3 < -0.5 * height) return false;
            return isInCircle(deltaX1, deltaX2, offset);
    }

    VF_LOG_CRITICAL("Unknown rotational axis in Cylinder.");
    return false;
}

void Cylinder::changeSizeByDelta(double delta)
{
    this->radius += delta;
    this->height += 2 * delta;
}
//! \}
