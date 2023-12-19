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
#ifndef OBJECT_H
#define OBJECT_H


#include "grid/Cell.h"
#include "global.h"

class GridImp;
struct Vertex;

class Object
{
public:
    virtual ~Object() = default;
    virtual SPtr<Object> clone() const = 0;

    virtual double getX1Centroid() const = 0;
    virtual double getX1Minimum() const  = 0;
    virtual double getX1Maximum() const  = 0;

    virtual double getX2Centroid() const = 0;
    virtual double getX2Minimum() const  = 0;
    virtual double getX2Maximum() const  = 0;

    virtual double getX3Centroid() const = 0;
    virtual double getX3Minimum() const  = 0;
    virtual double getX3Maximum() const  = 0;


    virtual void changeSizeByDelta(double delta) = 0;


    virtual bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) = 0;

    virtual bool isCellInObject(const Cell& cell) {
        for (const auto point : cell)
        {
            const bool isInObject = isPointInObject(point.x, point.y, point.z, 0.0, 0.0);
            if (!isInObject)
                return false;
        }
        return true;
    }

    virtual void findInnerNodes(SPtr<GridImp> grid);

    virtual int getIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnObject, real &qVal);
};


#endif

//! \}
