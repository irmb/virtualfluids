#ifndef LOG_FILE_DATA_ASSISTANT_STRATEGY_FACTORY_H
#define LOG_FILE_DATA_ASSISTANT_STRATEGY_FACTORY_H

#include "Simulation/BasicSimulation.h"

#include <memory>

class LogFileDataAssistantStrategy;

class LogFileDataAssistantStrategyFactory
{
public:
    virtual std::shared_ptr<LogFileDataAssistantStrategy> makeLogFileDataAssistantStrategy(BasicSimulation sim) = 0;
};
#endif 