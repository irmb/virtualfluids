#ifndef LOG_FILE_DATA_ASSISTANT_STRATEGY_H
#define LOG_FILE_DATA_ASSISTANT_STRATEGY_H

#include<memory>
#include <string>

class LogFileData;

class LogFileDataAssistantStrategy
{
public:
    virtual std::string getSimulationName() = 0;
    virtual bool checkSimulationParameter(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2) = 0;
};
#endif