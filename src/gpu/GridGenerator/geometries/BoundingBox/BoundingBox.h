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
#ifndef BoundingBox_h
#define BoundingBox_h

#include <vector>

#include "global.h"

struct Vertex;
struct Triangle;


class BoundingBox
{
public:
    real minX;
    real maxX;
    real minY;
    real maxY;
    real minZ;
    real maxZ;

    BoundingBox(real minX, real maxX, real minY, real maxY, real minZ, real maxZ);
    BoundingBox() = default;

public:
    static BoundingBox makeInvalidMinMaxBox();

    void setMinMax(const Triangle &t);
    void print() const;

    bool isInside(const Triangle &t) const;
    bool isInside(const real x, const real y, const real z) const;
    bool intersect(const Triangle &t) const;

    std::vector<std::vector<Vertex> > getIntersectionPoints(const BoundingBox &b) const;
    bool intersect(const BoundingBox &box) const;

    bool operator==(const BoundingBox &box) const;
    
    void extend(real delta);

private:
    bool isInside(const Vertex &v) const;
    void getPoints(Vertex v[8]) const;

};

#endif


//! \}
