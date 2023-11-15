#ifndef LOGFILE_DATA_ASSISTANT_H
#define LOGFILE_DATA_ASSISTANT_H

#include "Simulation/BasicSimulation.h"

#include <memory>
#include <vector>
#include <string>

enum DataCombination{ EqualSimulationsForDifferentKernels , EqualKernelSimulationsForDifferentViscosities};

class LogFileData;
class LogFileDataGroup;
class LogFileDataAssistantStrategy;

class LogFileDataAssistant
{
public:
    virtual std::vector<std::shared_ptr<LogFileDataGroup> > findDataCombination(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy, DataCombination combination) = 0;

};
#endif