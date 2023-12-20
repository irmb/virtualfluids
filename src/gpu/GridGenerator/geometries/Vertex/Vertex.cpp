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
#include "Vertex.h"

#include "utilities/math/Math.h"

Vertex::Vertex(real x, real y, real z) : x(x), y(y), z(z){}
Vertex::Vertex() { x = 0.0f; y = 0.0f; z = 0.0f; }

 real Vertex::getEuclideanDistanceTo(const Vertex &w) const
{
    return vf::Math::sqrtReal((x - w.x)*(x - w.x) + (y - w.y)*(y - w.y) + (z - w.z)*(z - w.z));
}

Vertex Vertex::operator-(const Vertex &v) const
{
    return Vertex(x - v.x, y - v.y, z - v.z);
}

Vertex Vertex::operator+(const Vertex &v) const
{
    return Vertex(this->x + v.x, this->y + v.y, this->z + v.z);
}

Vertex Vertex::operator*(const real& value) const
{
    return Vertex(value * this->x, value * this->y, value * this->z);
}

Vertex Vertex::operator/(const real& value) const
{ 
    return *this * ((real)1.0 / value); 
}

real Vertex::operator*(const Vertex &w) const
{
    return x*w.x + y*w.y + z*w.z;
}

struct Vertex Vertex::crossProduct(const Vertex &w) const
{
    real a = y*w.z - z*w.y;
    real b = z*w.x - x*w.z;
    real c = x*w.y - y*w.x;
    return Vertex(a, b, c);
}

real Vertex::length() const 
{
    return vf::Math::sqrtReal(x * x + y * y + z * z);
}

void Vertex::normalize()
{
    real len = length();

    if (len > EPSILON)
    {
        real invLen = 1.0f / len;
        x *= invLen;
        y *= invLen;
        z *= invLen;
    }
}

real Vertex::getMagnitude() const
{
    real temp = x*x + y*y + z*z;
    return vf::Math::sqrtReal(temp);
}

int Vertex::isEqual(const Vertex &w) const
{
    return vf::Math::equal(x, w.x) && vf::Math::equal(y, w.y) && vf::Math::equal(z, w.z);
}

real Vertex::getInnerAngle(const Vertex &w) const
{
    if (isEqual(w))
        return 0.0;
    if (this->getMagnitude() == 0 || w.getMagnitude() == 0)
        return 0.0;

    real mag = this->getMagnitude() * w.getMagnitude();
    real skal = *this * w;
    if (mag - std::abs(skal) < 0.0001)
        return 0.0f;
    return  vf::Math::acosReal(skal / mag) * 180.0f / vf::Math::acosReal(-1.0f); // acos(-1.0f) = PI 
}

void Vertex::print() const
{
    printf("(%2.8f,%2.8f,%2.8f)\n", x, y, z);
}

void Vertex::print(std::ostream &ost) const
{
    ost.write((char*)&x, 4);
    ost.write((char*)&y, 4);
    ost.write((char*)&z, 4);
}

void Vertex::printFormatted(std::ostream &ost) const
{
    ost << x << " " << y << " " << z;
}



bool Vertex::operator==(const Vertex &v) const
{
    return vf::Math::equal(x, v.x) && vf::Math::equal(y, v.y) && vf::Math::equal(z, v.z);
}


bool Vertex::isXbetween(real min, real max) const
{
    return x >= min && x <= max;
}

bool Vertex::isYbetween(real min, real max) const
{
    return y >= min && y <= max;
}

bool Vertex::isZbetween(real min, real max) const
{
    return z >= min && z <= max;
}

void Vertex::setMinMax(real & minX, real & maxX, real & minY, real & maxY, real & minZ, real & maxZ, const Vertex & v1, const Vertex & v2, const Vertex & v3)
{
    calculateMinMax(v1.x, v2.x, v3.x, minX, maxX);
    calculateMinMax(v1.y, v2.y, v3.y, minY, maxY);
    calculateMinMax(v1.z, v2.z, v3.z, minZ, maxZ);
}


real getMinimum(const real &value1, const real &value2)
{
    return value1 < value2 ? value1 : value2;
}

real getMaximum(const real &value1, const real &value2)
{
    return value1 > value2 ? value1 : value2;
}


void Vertex::calculateMinMax(const real &value1, const real &value2, const real &value3, real &min, real &max)
{
    
    real newMinimum = value1;
    newMinimum = getMinimum(value2, newMinimum);
    newMinimum = getMinimum(value3, newMinimum);

    real newMaximum = value1;
    newMaximum = getMaximum(value2, newMaximum);
    newMaximum = getMaximum(value3, newMaximum);

    min = newMinimum;
    max = newMaximum;
}



//! \}
