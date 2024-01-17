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
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "ShearWaveSimulationParameter.h"

#include "Simulations/ShearWave/InitialConditions/InitialConditionShearWave.h"
#include "Simulations/ShearWave/ShearWaveParameterStruct.h"

#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameter> ShearWaveSimulationParameter::getNewInstance(std::string kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
    return std::shared_ptr<SimulationParameter>(new ShearWaveSimulationParameter(kernel, viscosity, parameterStruct, gridInfo));
}

ShearWaveSimulationParameter::ShearWaveSimulationParameter(std::string kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> parameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernel, viscosity, parameterStruct->basicSimulationParameter, gridInfo)
{
    this->timeStepLength = parameterStruct->basicTimeStepLength * (gridInfo->lx / l0)*(gridInfo->lx / l0);

    if (parameterStruct->ux > parameterStruct->uz)
        this->maxVelocity = parameterStruct->ux / (lx / l0);
    else
        this->maxVelocity = parameterStruct->uz / (lx / l0);
 
    std::string kernelName = kernel;

    std::ostringstream oss;
    oss << parameterStruct->vtkFilePath << "/ShearWave/Viscosity_" << viscosity << "/ux_" << parameterStruct->ux << "_uz_" << parameterStruct->uz << "/" << kernelName << "/grid" << lx;
    generateFileDirectionInMyStystem(oss.str());
    this->filePath = oss.str();
}
//! \}
