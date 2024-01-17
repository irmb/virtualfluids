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
#include "SimulationParameterTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/InitialConditions/InitialConditionTaylorGreenVortexUz.h"
#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"

#include <sstream>

std::shared_ptr<SimulationParameterTaylorGreenUz> SimulationParameterTaylorGreenUz::getNewInstance(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
{
    return std::shared_ptr<SimulationParameterTaylorGreenUz>(new SimulationParameterTaylorGreenUz(kernel, viscosity, tgvParameterStruct, gridInfo));
}

SimulationParameterTaylorGreenUz::SimulationParameterTaylorGreenUz(std::string kernel, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> tgvParameterStruct, std::shared_ptr<GridInformationStruct> gridInfo)
:SimulationParameterImp(kernel, viscosity, tgvParameterStruct->basicSimulationParameter, gridInfo)
{
    this->timeStepLength = tgvParameterStruct->basicTimeStepLength * (gridInfo->lz / l0)*(gridInfo->lz / l0);
    this->maxVelocity = tgvParameterStruct->uz / (lz / l0);

    std::string kernelName = kernel;

    std::ostringstream oss;
    oss << tgvParameterStruct->vtkFilePath << "/TaylorGreenVortex Uz/Viscosity_" << viscosity << "/uz_" << tgvParameterStruct->uz << "_amplitude_" << tgvParameterStruct->amplitude << "/" << kernelName << "/grid" << lx;
    generateFileDirectionInMyStystem(oss.str());
    this->filePath = oss.str();
}
//! \}
