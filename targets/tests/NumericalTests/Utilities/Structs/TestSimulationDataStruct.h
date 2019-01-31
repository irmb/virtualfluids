#ifndef TEST_SIMULATION_DATA_STRUCT_H
#define TEST_SIMULATION_DATA_STRUCT_H

#include <memory>

class SimulationParameter;
class SimulationInfo;
class AnalyticalResults;

struct TestSimulationDataStruct
{
	std::shared_ptr< SimulationParameter> simParameter;
	std::shared_ptr< SimulationInfo> simInformation;
	std::shared_ptr< AnalyticalResults> analyticalResult;
};
#endif 