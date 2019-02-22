#ifndef LOGFILE_DATA_ASSISTANT_H
#define LOGFILE_DATA_ASSISTANT_H

#include "Simulation/BasicSimulation.h"

#include <memory>
#include <vector>

enum DataCombination{ EqualSimulationsForDifferentKernels , EqualKernelSimulationsForDifferentViscosities};

class LogFileData;
class LogFileDataGroup;

class LogFileDataAssistant
{
public:
	virtual std::vector<std::shared_ptr<LogFileDataGroup> > findDataCombination(std::vector<std::shared_ptr<LogFileData> > allLogFileData, BasicSimulation simulation, DataCombination combination) = 0;

};
#endif