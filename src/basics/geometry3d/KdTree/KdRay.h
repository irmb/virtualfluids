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
//! \author Soeren Textor, Sebastian Bindick
//=======================================================================================
#ifndef KDRAY_H
#define KDRAY_H

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbMath.h>

namespace Kd
{
/*
 * Ray class, for use with the optimized ray-box intersection test
 * described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 *
 */
template <typename T>
class Ray
{
public:
    Ray(const T &originX, const T &originY, const T &originZ, const T &directionX, const T &directionY,
        const T &directionZ)
    {
        this->originX = originX;
        this->originY = originY;
        this->originZ = originZ;

        // normierung (fuer ray-triangle-intersection)
        T oneOverLength =
            T(1.0 / std::sqrt(directionX * directionX + directionY * directionY + directionZ * directionZ));

        this->directionX = directionX * oneOverLength;
        this->directionY = directionY * oneOverLength;
        this->directionZ = directionZ * oneOverLength;

        this->inv_directionX = T(1.0 / this->directionX); // ACHTUNG: BEWUSST KEINE ==0 Abfrage
        this->inv_directionY = T(1.0 / this->directionY); // Alg verwendet exlitzit INF
        this->inv_directionZ = T(1.0 / this->directionZ);

        if (this->inv_directionX < 0.0)
            this->signX = 1;
        else
            this->signX = 0;
        if (this->inv_directionY < 0.0)
            this->signY = 1;
        else
            this->signY = 0;
        if (this->inv_directionZ < 0.0)
            this->signZ = 1;
        else
            this->signZ = 0;
    }

    T originX;
    T originY;
    T originZ;

    T directionX;
    T directionY;
    T directionZ;

    T inv_directionX;
    T inv_directionY;
    T inv_directionZ;

    int signX;
    int signY;
    int signZ;
};
} // namespace Kd

#endif // KDRAY_H

//! \}
