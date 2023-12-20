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
#include "SimulationFileNames.h"

const std::string fileEnding = ".dat";

const std::string simulationFileNames::coordX = "coordX" + fileEnding;
const std::string simulationFileNames::coordY = "coordY" + fileEnding;
const std::string simulationFileNames::coordZ = "coordZ" + fileEnding;
const std::string simulationFileNames::neighborX = "neighborX" + fileEnding;
const std::string simulationFileNames::neighborY = "neighborY" + fileEnding;
const std::string simulationFileNames::neighborZ = "neighborZ" + fileEnding;
const std::string simulationFileNames::neighborWSB = "neighborWSB" + fileEnding;
const std::string simulationFileNames::geoVec = "geoVec" + fileEnding;

const std::string simulationFileNames::scaleCFC = "scaleCFC" + fileEnding;
const std::string simulationFileNames::scaleCFF = "scaleCFF" + fileEnding;
const std::string simulationFileNames::scaleFCC = "scaleFCC" + fileEnding;
const std::string simulationFileNames::scaleFCF = "scaleFCF" + fileEnding;

const std::string simulationFileNames::offsetVecCF = "offsetVecCF" + fileEnding;
const std::string simulationFileNames::offsetVecFC = "offsetVecFC" + fileEnding;
        
const std::string simulationFileNames::geomBoundaryQ = "geomBoundaryQs" + fileEnding;
const std::string simulationFileNames::geomBoundaryValues = "geomBoundaryValues" + fileEnding;
            
const std::string simulationFileNames::topBoundaryQ = "topBoundaryQs" + fileEnding;
const std::string simulationFileNames::topBoundaryValues = "topBoundaryValues" + fileEnding;
                
const std::string simulationFileNames::bottomBoundaryQ = "bottomBoundaryQs" + fileEnding;
const std::string simulationFileNames::bottomBoundaryValues = "bottomBoundaryValues" + fileEnding;
            
const std::string simulationFileNames::frontBoundaryQ = "frontBoundaryQs" + fileEnding;
const std::string simulationFileNames::frontBoundaryValues = "frontBoundaryValues" + fileEnding;
                 
const std::string simulationFileNames::backBoundaryQ = "backBoundaryQs" + fileEnding;
const std::string simulationFileNames::backBoundaryValues = "backBoundaryValues" + fileEnding;
                  
const std::string simulationFileNames::inletBoundaryQ = "inletBoundaryQs" + fileEnding;
const std::string simulationFileNames::inletBoundaryValues = "inletBoundaryValues" + fileEnding;
                
const std::string simulationFileNames::outletBoundaryQ = "outletBoundaryQs" + fileEnding;
const std::string simulationFileNames::outletBoundaryValues = "outletBoundaryValues" + fileEnding;

const std::string simulationFileNames::numberNodes = "numberNodes" + fileEnding;
const std::string simulationFileNames::LBMvsSI = "LBMvsSI" + fileEnding;


//! \}
