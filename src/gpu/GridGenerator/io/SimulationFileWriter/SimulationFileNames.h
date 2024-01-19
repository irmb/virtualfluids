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
//! \addtogroup gpu_io io
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef SimulationFileNames_H
#define SimulationFileNames_H

#include <string>

#include "global.h"

struct simulationFileNames
{
    static const std::string coordX;
    static const std::string coordY;
    static const std::string coordZ;
    static const std::string neighborX;
    static const std::string neighborY;
    static const std::string neighborZ;
    static const std::string neighborWSB;
    static const std::string geoVec;

    static const std::string scaleCFC;
    static const std::string scaleCFF;
    static const std::string scaleFCC;
    static const std::string scaleFCF;

    static const std::string offsetVecCF;
    static const std::string offsetVecFC;
    
    static const std::string geomBoundaryQ;
    static const std::string geomBoundaryValues;
    
    static const std::string topBoundaryQ;
    static const std::string topBoundaryValues;
    
    static const std::string bottomBoundaryQ;
    static const std::string bottomBoundaryValues;
    
    static const std::string frontBoundaryQ;
    static const std::string frontBoundaryValues;
    
    static const std::string backBoundaryQ;
    static const std::string backBoundaryValues;
    
    static const std::string inletBoundaryQ;
    static const std::string inletBoundaryValues;
    
    static const std::string outletBoundaryQ;
    static const std::string outletBoundaryValues;

    static const std::string numberNodes;
    static const std::string LBMvsSI;
};


#endif

//! \}
