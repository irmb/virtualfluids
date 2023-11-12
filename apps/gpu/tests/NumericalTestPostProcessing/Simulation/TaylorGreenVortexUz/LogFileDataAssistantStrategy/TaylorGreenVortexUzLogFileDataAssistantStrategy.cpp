#include "TaylorGreenVortexUzLogFileDataAssistantStrategy.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Simulation/TaylorGreenVortexUz/LogFileData/TaylorGreenVortexUzLogFileData.h"

std::shared_ptr<LogFileDataAssistantStrategy> TaylorGreenVortexUzLogFileDataAssistantStrategy::getNewInstance()
{
    return std::shared_ptr<LogFileDataAssistantStrategy>(new TaylorGreenVortexUzLogFileDataAssistantStrategy());
}

std::string TaylorGreenVortexUzLogFileDataAssistantStrategy::getSimulationName()
{
    return simName;
}

bool TaylorGreenVortexUzLogFileDataAssistantStrategy::checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
    if (!equalDouble(logFileData1->getTaylorGreenVortexUzLogFileData()->getUz().at(0), logFileData2->getTaylorGreenVortexUzLogFileData()->getUz().at(0)))
        return false;
    if (!equalDouble(logFileData1->getTaylorGreenVortexUzLogFileData()->getAmplitude().at(0), logFileData2->getTaylorGreenVortexUzLogFileData()->getAmplitude().at(0)))
        return false;
    if (logFileData1->getTaylorGreenVortexUzLogFileData()->getL0().at(0) != logFileData2->getTaylorGreenVortexUzLogFileData()->getL0().at(0))
        return false;

    return true;
}

TaylorGreenVortexUzLogFileDataAssistantStrategy::TaylorGreenVortexUzLogFileDataAssistantStrategy()
{
    this->simName = "TaylorGreenVortexUz";
}
