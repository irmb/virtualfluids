#ifndef SIMULATION_DATA_STRUCT_H
#define SIMULATION_DATA_STRUCT_H

#include "Utilities/Structs/TestSimulationDataStruct.h"

#include <vector>
#include <memory>

struct SimulationDataStruct
{
	std::vector<std::shared_ptr<TestSimulationDataStruct> > testSimData;

	std::shared_ptr<SimulationLogFileInformation> logFileInformation;
	bool simGroupRun;
};
#endif 