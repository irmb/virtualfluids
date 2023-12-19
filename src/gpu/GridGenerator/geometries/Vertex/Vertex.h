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
#ifndef VERTEX_H
#define VERTEX_H

#include <stdio.h>
#include <memory>
#include <ostream>

#include "gpu/GridGenerator/global.h"

struct Vertex
{
public:
    real x, y, z;

    Vertex(real x, real y, real z);
    Vertex();

    real getEuclideanDistanceTo(const Vertex &w) const;
    Vertex operator-(const Vertex &v) const;
    Vertex operator+(const Vertex &v) const;
    Vertex operator*(const real& value) const;
    Vertex operator/(const real& value) const;

    real operator*(const Vertex &w) const;
    struct Vertex crossProduct(const Vertex &w) const;
    real length() const;
    void normalize();
    real getMagnitude() const;
    int isEqual(const Vertex &w) const;
    real getInnerAngle(const Vertex &w) const;

    bool operator==(const Vertex &v) const;

    bool isXbetween(real min, real max) const;
    bool isYbetween(real min, real max) const;
    bool isZbetween(real min, real max) const;

    static void setMinMax(real &minX, real &maxX, real &minY, real &maxY, real &minZ, real &maxZ, const Vertex &v1, const Vertex &v2, const Vertex &v3); 
    static void calculateMinMax(const real &value1, const real &value2, const real &value3, real &min, real &max);

    void print() const;
    void print(std::ostream &ost) const;
    void printFormatted(std::ostream &ost) const;

};



#endif

//! \}
