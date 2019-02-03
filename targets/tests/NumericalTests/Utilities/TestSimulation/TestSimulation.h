#ifndef TEST_SIMULATION_H
#define TEST_SIMULATION_H

#include "Utilities\NumericalTestSuite\NumericalTestSuite.h"

#include <memory>
#include <string>

class AnalyticalResults;
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
	//beachte Interface Segregation
	//TimeStopper

	virtual std::shared_ptr<SimulationParameter> getSimulationParameter() = 0;
	virtual std::shared_ptr<DataWriter> getDataWriter() = 0;
	virtual std::shared_ptr<InitialCondition> getInitialCondition() = 0;
	virtual std::shared_ptr<TimeTracking> getTimeTracking() = 0;
	virtual void setParameter(std::shared_ptr<Parameter> para) = 0;

	virtual bool getSimulationRun() = 0;

	virtual void makeSimulationHeadOutput() = 0;
	virtual void startPostProcessing() = 0;
};
#endif