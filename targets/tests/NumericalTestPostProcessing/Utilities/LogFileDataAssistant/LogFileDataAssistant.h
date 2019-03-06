#ifndef LOGFILE_DATA_ASSISTANT_H
#define LOGFILE_DATA_ASSISTANT_H

#include "Simulation/BasicSimulation.h"

#include <memory>
#include <vector>

class LogFileData;
class LogFileDataGroup;

class LogFileDataAssistant
{
public:
	virtual std::vector<std::shared_ptr<LogFileDataGroup> > findEqualSimulationsForDifferentKernels(std::vector<std::shared_ptr<LogFileData> > allLogFileData, BasicSimulation simulation) = 0;
	virtual std::vector<std::shared_ptr<LogFileDataGroup> > findEqualKernelSimulationsForDifferentViscosities(std::vector<std::shared_ptr<LogFileData> > allLogFileData, BasicSimulation simulation) = 0;
};
#endif