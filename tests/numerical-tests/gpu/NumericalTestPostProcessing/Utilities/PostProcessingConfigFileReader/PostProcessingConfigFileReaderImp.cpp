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
#include "PostProcessingConfigFileReaderImp.h"

#include <basics/config/ConfigurationFile.h>
#include "StringUtilities/StringUtil.h"

#include "Utilities/PostProcessingConfigData/PostProcessingConfigDataImp.h"

#include <fstream>


std::shared_ptr<PostProcessingConfigFileReader> PostProcessingConfigFileReaderImp::getNewInstance()
{
    return std::shared_ptr<PostProcessingConfigFileReader>(new PostProcessingConfigFileReaderImp());
}

std::shared_ptr<PostProcessingConfigData> PostProcessingConfigFileReaderImp::readConfigFile(std::string filePath)
{
    auto input = std::make_shared<vf::basics::ConfigurationFile>();
    input->load(filePath);

    std::vector<BasicSimulation> simulation;
    std::vector<Assistant> assistants;
    std::vector<DataCombination> combination;

    if(StringUtil::toBool(input->getValue<std::string>("ShearWave")))
        simulation.push_back(ShearWave);

    if (StringUtil::toBool(input->getValue<std::string>("TaylorGreenVortexUx")))
        simulation.push_back(TaylorGreenVortexUx);

    if (StringUtil::toBool(input->getValue<std::string>("TaylorGreenVortexUz")))
        simulation.push_back(TaylorGreenVortexUz);

    if (StringUtil::toBool(input->getValue<std::string>("Phi")))
        assistants.push_back(Phi);

    if (StringUtil::toBool(input->getValue<std::string>("Ny")))
        assistants.push_back(Ny);

    if (StringUtil::toBool(input->getValue<std::string>("L2Norm")))
        assistants.push_back(L2Norm);

    if (StringUtil::toBool(input->getValue<std::string>("L2Norm_BetweenKernels")))
        assistants.push_back(L2NormBetweenKernels);

    if (StringUtil::toBool(input->getValue<std::string>("TimeOutput")))
        assistants.push_back(Time);


    if (StringUtil::toBool(input->getValue<std::string>("EqualSimulationsForDifferentKernels")))
        combination.push_back(EqualSimulationsForDifferentKernels);

    if (StringUtil::toBool(input->getValue<std::string>("EqualKernelSimulationsForDifferentViscosities")))
        combination.push_back(EqualKernelSimulationsForDifferentViscosities);

    std::shared_ptr<PostProcessingConfigDataImp> data = PostProcessingConfigDataImp::getNewInstance();

    data->setAssistants(assistants);
    data->setSimulations(simulation);
    data->setDataCombinations(combination);

    data->setLogFilesPath(input->getValue<std::string>("LogFilesPath"));
    data->setMathematicaFilePath(input->getValue<std::string>("MathematicaFilePath"));
    
    return data;
}

//! \}
