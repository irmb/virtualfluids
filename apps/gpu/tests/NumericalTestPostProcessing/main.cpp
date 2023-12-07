#include "Simulation/BasicSimulation.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"
#include "Utilities/LogFileReader/LogFileReader.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantImp.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantStrategy/LogFileDataAssistantStrategyFactory/LogFileDataAssistantStrategyFactoryImp.h"
#include "Utilities/MathematicaFile/MathematicaFile.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactoryImp.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistantFactory/MathematicaAssistantFactoryImp.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistant.h"

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
