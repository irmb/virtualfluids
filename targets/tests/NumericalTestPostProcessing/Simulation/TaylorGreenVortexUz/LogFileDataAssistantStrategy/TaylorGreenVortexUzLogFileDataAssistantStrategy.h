#ifndef LOG_FILE_DATA_ASSISTANT_STRATEGY_TGVUZ_H
#define LOG_FILE_DATA_ASSISTANT_STRATEGY_TGVUZ_H

#include "Utilities/LogFileDataAssistant/LogFileDataAssistantStrategy/LogFileDataAssistantStrategyImp.h"

class TaylorGreenVortexUzLogFileDataAssistantStrategy : public LogFileDataAssistantStrategyImp
{
public:
	static std::shared_ptr<LogFileDataAssistantStrategy> getNewInstance();

	std::string getSimulationName();
	bool checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);

private:
	TaylorGreenVortexUzLogFileDataAssistantStrategy();

	std::string simName;
};
#endif