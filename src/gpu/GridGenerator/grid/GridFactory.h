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
#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H

#include "global.h"

#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "grid/GridImp.h"

namespace vf::gpu {

enum class TriangularMeshDiscretizationMethod
{
    RAYCASTING, POINT_IN_OBJECT, POINT_UNDER_TRIANGLE
#if defined(VF_HAS_FAST_WINDING)
    , FAST_WINDING
#endif
};

class GridFactory
{
public:
    static SPtr<GridFactory> make()
    {
        return std::make_shared<GridFactory>();
    }

    SPtr<Grid> makeGrid(SPtr<Object> gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level, const std::string& d3Qxx = "D3Q27")
    {
        SPtr<GridImp> grid;
        
        grid = GridImp::makeShared(gridShape, startX, startY, startZ, endX, endY, endZ, delta, d3Qxx, level);

        if (triangularMeshDiscretizationStrategy)
            grid->setTriangularMeshDiscretizationStrategy(triangularMeshDiscretizationStrategy);
        else
            grid->setTriangularMeshDiscretizationStrategy(std::make_shared<PointInObjectDiscretizationStrategy>());

        return grid;
    }

    void setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod triangularMeshDiscretizationMethod)
    {
        switch (triangularMeshDiscretizationMethod)
        {
        case TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE:
            triangularMeshDiscretizationStrategy = std::make_shared<PointUnderTriangleStrategy>();
            break;
        case TriangularMeshDiscretizationMethod::RAYCASTING:
            triangularMeshDiscretizationStrategy = std::make_shared<RayCastingDiscretizationStrategy>();
            break;
        case TriangularMeshDiscretizationMethod::POINT_IN_OBJECT:
            triangularMeshDiscretizationStrategy = std::make_shared<PointInObjectDiscretizationStrategy>();
            break;
#if defined(VF_HAS_FAST_WINDING)
        case TriangularMeshDiscretizationMethod::FAST_WINDING:
            triangularMeshDiscretizationStrategy = std::make_shared<FastWindingDiscretizationStrategy>();
            break;
#endif
        }
    }

    void setTriangularMeshDiscretizationStrategy(
        SPtr<TriangularMeshDiscretizationStrategy> strategy)
    {
        triangularMeshDiscretizationStrategy = strategy;
    }

private:
    SPtr<TriangularMeshDiscretizationStrategy> triangularMeshDiscretizationStrategy;
};

}

#endif

//! \}
