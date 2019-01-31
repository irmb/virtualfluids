#ifndef TEST_SIMULATION_H
#define TEST_SIMULATION_H

#include "SimulationTime\SimulationTime.h"

#include <memory>

class AnalyticalResults;
class SimulationParameter;
class SimulationResults;
class SimulationObserver;
class DataWriter;
class Parameter;
class SimulationInfo;

class TestSimulation
{
public:
	//beachte Interface Segregation
	//TimeStopper

	virtual std::shared_ptr<SimulationParameter> getSimulationParameter() = 0;
	virtual std::shared_ptr<DataWriter> getDataWriter() = 0;

	virtual bool getSimulationRun() = 0;
	virtual std::shared_ptr<SimulationResults> getSimulationResults() = 0;
	virtual std::string getRunTimeOutput() = 0;


	virtual void registerSimulationObserver(std::shared_ptr< SimulationObserver> simObserver) = 0;
	virtual void makeSimulationHeadOutput() = 0;
	virtual void setSimulationStartTime() = 0;
	virtual void setSimulationEndTimeAndNotifyObserver() = 0;

	virtual void setParameter(std::shared_ptr<Parameter> para) = 0;

};
#endif