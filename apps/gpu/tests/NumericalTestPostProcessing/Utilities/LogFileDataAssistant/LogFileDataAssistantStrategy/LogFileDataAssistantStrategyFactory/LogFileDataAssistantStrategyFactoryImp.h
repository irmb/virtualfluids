#ifndef LOG_FILE_DATA_ASSISTANT_STRATEGY_FACTORY_IMP_H
#define LOG_FILE_DATA_ASSISTANT_STRATEGY_FACTORY_IMP_H

#include "LogFileDataAssistantStrategyFactory.h"

class LogFileDataAssistantStrategyFactoryImp : public LogFileDataAssistantStrategyFactory
{
public:
    static std::shared_ptr<LogFileDataAssistantStrategyFactory> getNewInstance();

    std::shared_ptr<LogFileDataAssistantStrategy> makeLogFileDataAssistantStrategy(BasicSimulation sim);
    
private:
    LogFileDataAssistantStrategyFactoryImp();
};
#endif 