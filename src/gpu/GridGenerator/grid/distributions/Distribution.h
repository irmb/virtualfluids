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
#ifndef Distribution_H
#define Distribution_H

#include <vector>
#include <string>

#include "gpu/GridGenerator/global.h"

#define DIR_END_MAX 27


struct Direction
{
    Direction()
    {
        dir[0] = 0;
        dir[1] = 0;
        dir[2] = 0;
    }

    Direction(int dx, int dy, int dz)
    {
        dir[0] = dx;
        dir[1] = dy;
        dir[2] = dz;
    }

    int operator[](uint dir) const
    {
        if (dir < 3)
            return this->dir[dir];
        return -99;
    }
private:
    int dir[3];
};

struct Distribution
{
    typedef Direction* iterator;
    typedef const Direction* const_iterator;

    real* f;
    std::vector<int> dirs;
    std::vector<Direction> directions;
    int dir_start;
    int dir_end;
    const char* name;
    unsigned long fSize;

    void setSize(uint size) {
        fSize = size * (dir_end + 1);
    }

    iterator begin() { return &directions[0]; }
    const_iterator begin() const { return &directions[0]; }
    iterator end() { return &directions[dir_end + 1]; }
    const_iterator end() const { return &directions[dir_end + 1]; }
};

class Grid;

class DistributionHelper
{
public:
    static Distribution getDistribution27();

    static Distribution getDistribution(std::string name);

public:
    static std::vector<std::vector<real> > getQsWithoutRowsWithOnlyZeroValues(const Grid &grid, const Distribution &d);
    static std::vector<std::vector<real> > getAllQsOnFluidNodes(const Grid &grid, const Distribution &d);
    static int getNeighborNodeIndexInGivenDirection(const Distribution &d, const Grid &grid, const int node, const int dir_index);
    static std::vector<std::vector<real> > getVectorWithoutRowsWithOnlyZeroValues(std::vector<std::vector<real> > qs);
    static void printQs(std::vector<std::vector<real> > qs, int decimalPlaces);
};

#endif

//! \}
