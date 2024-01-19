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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef CELL_H
#define CELL_H

#include "gpu/GridGenerator/global.h"

struct Point
{
    Point() : x(0.0), y(0.0), z(0.0) {}
    Point(real x, real y, real z) : x(x), y(y), z(z) {}
    real x, y, z;
};

class Cell
{
public:
    typedef Point* iterator;
    typedef const Point* const_iterator;

    Cell(real startX, real startY, real startZ, real delta)
    {
        points = new Point[size];
        points[0] = Point(startX, startY, startZ); // 0,0,0
        points[1] = Point(startX + delta, startY, startZ); // 1,0,0
        points[2] = Point(startX, startY + delta, startZ); // 0,1,0
        points[3] = Point(startX + delta, startY + delta, startZ); // 1,1,0

        points[4] = Point(startX, startY, startZ + delta); // 0,0,1
        points[5] = Point(startX + delta, startY, startZ + delta); // 1,0,1
        points[6] = Point(startX, startY + delta, startZ + delta); // 0,1,1
        points[7] = Point(startX + delta, startY + delta, startZ + delta); // 1,1,1
    }

    ~Cell()
    {
        delete[] points;
    }

    iterator begin() { return &points[0]; }
    const_iterator begin() const { return &points[0]; }
    iterator end() { return &points[size]; }
    const_iterator end() const { return &points[size]; }

private:
    Point* points;
    uint size = 8;

};



#endif

//! \}
