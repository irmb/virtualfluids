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
#include "Cuboid.h"

#include "PointerDefinitions.h"
#include "utilities/math/Math.h"

Cuboid::Cuboid(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b)
    : minX1(x1a), minX2(x2a), minX3(x3a), maxX1(x1b), maxX2(x2b), maxX3(x3b)
{

}

SPtr<Object> Cuboid::clone() const
{
    return std::make_shared<Cuboid>(minX1, minX2, minX3, maxX1, maxX2, maxX3);
}

double Cuboid::getX1Centroid() const
{
   return getCenter(minX1, maxX1);
}

double Cuboid::getX1Minimum() const
{
    return getMinimum(minX1, maxX1);
}

double Cuboid::getX1Maximum() const
{
    return getMaximum(minX1, maxX1);
}

double Cuboid::getX2Centroid() const
{
    return getCenter(minX2, maxX2);
}

double Cuboid::getX2Minimum() const
{
    return getMinimum(minX2, maxX2);
}    

double Cuboid::getX2Maximum() const
{
    return getMaximum(minX2, maxX2);
}

double Cuboid::getX3Centroid() const
{
    return getCenter(minX3, maxX3);
}

double Cuboid::getX3Minimum() const
{
    return getMinimum(minX3, maxX3);
}    

double Cuboid::getX3Maximum() const
{
    return getMaximum(minX3, maxX3);
}

double Cuboid::getCenter(double x1, double x2)
{
    return 0.5 * (x1 + x2);
}

double Cuboid::getMinimum(double x1, double x2)
{
    return (x1 < x2 ? x1 : x2);
}

double Cuboid::getMaximum(double x1, double x2)
{
    return (x1 > x2 ? x1 : x2);
}

bool Cuboid::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
    //false, if 'not in Object' or 'on Boundary'!
    if (vf::Math::lessEqual(   (real)x1, (real)this->getX1Minimum() + (real)minOffset))  return false;
    if (vf::Math::lessEqual(   (real)x2, (real)this->getX2Minimum() + (real)minOffset))  return false;
    if (vf::Math::lessEqual(   (real)x3, (real)this->getX3Minimum() + (real)minOffset))  return false;
    if (vf::Math::greaterEqual((real)x1, (real)this->getX1Maximum() - (real)maxOffset))  return false;
    if (vf::Math::greaterEqual((real)x2, (real)this->getX2Maximum() - (real)maxOffset))  return false;
    if (vf::Math::greaterEqual((real)x3, (real)this->getX3Maximum() - (real)maxOffset))  return false;

    return true;
}


bool Cuboid::isOn(const real& coord, const real& plane1, const real& plane2)
{
    return  vf::Math::equal(coord, plane1) || vf::Math::equal(coord, plane2);
}

bool Cuboid::isBetween(const real& coord, const real& start, const real& end)
{
    return  vf::Math::greaterEqual(coord, start) && vf::Math::lessEqual(coord, end);
}


void Cuboid::changeSizeByDelta(double delta)
{
    this->minX1 -= delta;
    this->minX2 -= delta;
    this->minX3 -= delta;
                   
    this->maxX1 += delta;
    this->maxX2 += delta;
    this->maxX3 += delta;
}

//! \}
