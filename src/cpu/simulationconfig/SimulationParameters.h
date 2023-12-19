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
#ifndef VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H
#define VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H

#include <array>
#include <memory>
#include <geometry3d/GbPoint3D.h>

struct PhysicalParameters
{
    double latticeViscosity{};
    double bulkViscosityFactor{1};
};

struct BoundingBox
{
    const double minX1;
    const double minX2;
    const double minX3;
    const double maxX1;
    const double maxX2;
    const double maxX3;

    BoundingBox(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3) :
            minX1(minX1),
            minX2(minX2),
            minX3(minX3),
            maxX1(maxX1),
            maxX2(maxX2),
            maxX3(maxX3)
    {}
};

struct GridParameters
{
    std::array<int, 3> numberOfNodesPerDirection{1, 1, 1};
    std::array<int, 3> blocksPerDirection{1, 1, 1};
    int referenceDirectionIndex{};
    double nodeDistance{1};
    bool periodicBoundaryInX1{};
    bool periodicBoundaryInX2{};
    bool periodicBoundaryInX3{};

    std::shared_ptr<BoundingBox> boundingBox()
    {
        return std::make_shared<BoundingBox>(
                0, 0, 0,
                numberOfNodesPerDirection[0] * nodeDistance,
                numberOfNodesPerDirection[1] * nodeDistance,
                numberOfNodesPerDirection[2] * nodeDistance
        );
    }
};

struct RuntimeParameters
{
    int numberOfTimeSteps{};
    int timeStepLogInterval{};
    int numberOfThreads{};
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H
