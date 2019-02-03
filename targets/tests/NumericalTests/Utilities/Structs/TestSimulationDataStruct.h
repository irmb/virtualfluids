#ifndef TEST_SIMULATION_DATA_STRUCT_H
#define TEST_SIMULATION_DATA_STRUCT_H

#include <memory>

class AnalyticalResults;
class InitialCondition;
class SimulationInfoImp;
class SimulationParameter;
class SimulationResults;
class TimeImp;

struct TestSimulationDataStruct
{
	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::shared_ptr<InitialCondition> initialCondition;
	std::shared_ptr<SimulationInfo> simInformation;
	std::shared_ptr<SimulationParameter> simParameter;
};
#endif 