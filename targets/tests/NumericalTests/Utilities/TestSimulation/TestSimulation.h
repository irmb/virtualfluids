#ifndef TEST_SIMULATION_H
#define TEST_SIMULATION_H

#include "Utilities\NumericalTestSuite\NumericalTestSuite.h"

#include <memory>
#include <string>

class InitialCondition;
class DataWriter;
class SimulationParameter;
class SimulationResults;
class SimulationObserver;
class TimeTracking;
class Parameter;

class TestSimulation : public NumericalTestSuite
{
public:
	virtual std::shared_ptr<SimulationParameter> getSimulationParameter() = 0;
	virtual std::shared_ptr<DataWriter> getDataWriter() = 0;
	virtual std::shared_ptr<InitialCondition> getInitialCondition() = 0;
	virtual std::shared_ptr<TimeTracking> getTimeTracking() = 0;
	virtual void setParameter(std::shared_ptr<Parameter> para) = 0;
};
#endif