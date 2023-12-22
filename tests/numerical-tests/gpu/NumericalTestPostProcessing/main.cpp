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
#include "Simulation/BasicSimulation.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantImp.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantStrategy/LogFileDataAssistantStrategyFactory/LogFileDataAssistantStrategyFactoryImp.h"
#include "Utilities/LogFileReader/LogFileReader.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistant.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistantFactory/MathematicaAssistantFactoryImp.h"
#include "Utilities/MathematicaFile/MathematicaFile.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactoryImp.h"

#include "Utilities/PostProcessingConfigFileReader/PostProcessingConfigFileReaderImp.h"
#include "Utilities/PostProcessingConfigData/PostProcessingConfigData.h"

#include "basics/DataTypes.h"

#include <memory>
#include <cfloat>
#include <cstdio>

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
    std::shared_ptr<PostProcessingConfigFileReader> reader = PostProcessingConfigFileReaderImp::getNewInstance();
    std::shared_ptr<PostProcessingConfigData> configData =  reader->readConfigFile(argv[1]);

    std::shared_ptr<LogFileReader> logFileReader = LogFileReader::getInstance();
    std::vector<std::shared_ptr<LogFileData> > logFileDataVector = logFileReader->readLogFilesInDirectoryToLogFileData(configData->getLogFilesPath());

    std::shared_ptr<MathematicaFile> aMathmaticaFile = MathematicaFile::getNewInstance(configData->getMathematicaFilePath());

    std::shared_ptr<LogFileDataAssistant> assistentLogFile = LogFileDataAssistantImp::getNewInstance();

    std::shared_ptr<LogFileDataAssistantStrategyFactory> assistentStrategyFactory = LogFileDataAssistantStrategyFactoryImp::getNewInstance();
    

    std::shared_ptr<MathematicaFunctionFactory> functionFactory = MathematicaFunctionFactoryImp::getNewInstance();
    std::shared_ptr<MathematicaAssistantFactory> assistantFactory = MathematicaAssistantFactoryImp::getNewInstance();
    std::vector<std::shared_ptr<MathematicaAssistant> > mathematicaAssistants = assistantFactory->makeMathematicaAssistants(configData->getAssistants(), functionFactory);

    for (uint sim = 0; sim < configData->getSimulations().size(); sim++) {
        for (uint comb = 0; comb < configData->getDataCombinations().size(); comb++) {
            std::shared_ptr<LogFileDataAssistantStrategy> strategy = assistentStrategyFactory->makeLogFileDataAssistantStrategy(configData->getSimulations().at(sim));
            std::vector<std::shared_ptr<LogFileDataGroup> > logFileDataSorted = assistentLogFile->findDataCombination(logFileDataVector, strategy, configData->getDataCombinations().at(comb));
            for (uint i = 0; i < logFileDataSorted.size(); i++) {
                for (uint j = 0; j < mathematicaAssistants.size(); j++)
                    mathematicaAssistants.at(j)->makeMathematicaOutput(logFileDataSorted.at(i), aMathmaticaFile);
            }
        }
    }
    
        
    aMathmaticaFile->finishFile();
    return 0;
}

//! \}
