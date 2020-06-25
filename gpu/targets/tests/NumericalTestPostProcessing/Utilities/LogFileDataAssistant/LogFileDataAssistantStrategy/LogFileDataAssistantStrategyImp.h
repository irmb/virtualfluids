#ifndef LOG_FILE_DATA_ASSISTANT_STRATEGY_IMP_H
#define LOG_FILE_DATA_ASSISTANT_STRATEGY_IMP_H

#include "LogFileDataAssistantStrategy.h"

class LogFileDataAssistantStrategyImp : public LogFileDataAssistantStrategy
{
public:
	virtual std::string getSimulationName() = 0;
	virtual bool checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2) = 0;

protected:
	bool equalDouble(double num1, double num2);
};
#endif