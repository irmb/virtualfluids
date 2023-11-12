#include "ShearWaveLogFileDataAssistantStrategy.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Simulation/ShearWave/LogFileData/ShearWaveLogFileData.h"

std::shared_ptr<LogFileDataAssistantStrategy> ShearWaveLogFileDataAssistantStrategy::getNewInstance()
{
    return std::shared_ptr<LogFileDataAssistantStrategy>(new ShearWaveLogFileDataAssistantStrategy());
}

std::string ShearWaveLogFileDataAssistantStrategy::getSimulationName()
{
    return simName;
}

bool ShearWaveLogFileDataAssistantStrategy::checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
    if (!equalDouble(logFileData1->getShearWaveLogFileData()->getUx().at(0), logFileData2->getShearWaveLogFileData()->getUx().at(0)))
        return false;
    if (!equalDouble(logFileData1->getShearWaveLogFileData()->getUz().at(0), logFileData2->getShearWaveLogFileData()->getUz().at(0)))
        return false;
    if (logFileData1->getShearWaveLogFileData()->getL0().at(0) != logFileData2->getShearWaveLogFileData()->getL0().at(0))
        return false;

    return true;
}

ShearWaveLogFileDataAssistantStrategy::ShearWaveLogFileDataAssistantStrategy()
{
    this->simName = "ShearWave";
}
