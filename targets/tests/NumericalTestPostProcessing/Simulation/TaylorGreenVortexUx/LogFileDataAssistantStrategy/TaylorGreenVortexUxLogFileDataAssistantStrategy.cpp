#include "TaylorGreenVortexUxLogFileDataAssistantStrategy.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Simulation/TaylorGreenVortexUx/LogFileData/TaylorGreenVortexUxLogFileData.h"

std::shared_ptr<LogFileDataAssistantStrategy> TaylorGreenVortexUxLogFileDataAssistantStrategy::getNewInstance()
{
	return std::shared_ptr<LogFileDataAssistantStrategy>(new TaylorGreenVortexUxLogFileDataAssistantStrategy());
}

std::string TaylorGreenVortexUxLogFileDataAssistantStrategy::getSimulationName()
{
	return simName;
}

bool TaylorGreenVortexUxLogFileDataAssistantStrategy::checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
	if (!equalDouble(logFileData1->getTaylorGreenVortexUxLogFileData()->getUx().at(0), logFileData2->getTaylorGreenVortexUxLogFileData()->getUx().at(0)))
		return false;
	if (!equalDouble(logFileData1->getTaylorGreenVortexUxLogFileData()->getAmplitude().at(0), logFileData2->getTaylorGreenVortexUxLogFileData()->getAmplitude().at(0)))
		return false;
	if (logFileData1->getTaylorGreenVortexUxLogFileData()->getL0().at(0) != logFileData2->getTaylorGreenVortexUxLogFileData()->getL0().at(0))
		return false;

	return true;
}

TaylorGreenVortexUxLogFileDataAssistantStrategy::TaylorGreenVortexUxLogFileDataAssistantStrategy()
{
	this->simName = "TaylorGreenVortexUx";
}
