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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "PostProcessingConfigDataImp.h"

PostProcessingConfigDataImp::PostProcessingConfigDataImp()
{
}

std::shared_ptr<PostProcessingConfigDataImp> PostProcessingConfigDataImp::getNewInstance()
{
    return std::shared_ptr<PostProcessingConfigDataImp>(new PostProcessingConfigDataImp());
}

std::vector<BasicSimulation> PostProcessingConfigDataImp::getSimulations()
{
    return simulations;
}

std::vector<Assistant> PostProcessingConfigDataImp::getAssistants()
{
    return assistants;
}

std::vector<DataCombination> PostProcessingConfigDataImp::getDataCombinations()
{
    return dataCombinations;
}

std::string PostProcessingConfigDataImp::getMathematicaFilePath()
{
    return mathematicaFilePath;
}

std::string PostProcessingConfigDataImp::getLogFilesPath()
{
    return logFilesPath;
}

void PostProcessingConfigDataImp::setSimulations(std::vector<BasicSimulation> simulations)
{
    this->simulations = simulations;
}

void PostProcessingConfigDataImp::setAssistants(std::vector<Assistant> assis)
{
    this->assistants = assis;
}

void PostProcessingConfigDataImp::setDataCombinations(std::vector<DataCombination> dataComb)
{
    this->dataCombinations = dataComb;
}

void PostProcessingConfigDataImp::setMathematicaFilePath(std::string mathematicaFilePath)
{
    this->mathematicaFilePath = mathematicaFilePath;
}

void PostProcessingConfigDataImp::setLogFilesPath(std::string logFilesPath)
{
    this->logFilesPath = logFilesPath;
}

//! \}
