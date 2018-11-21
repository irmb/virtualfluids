#ifndef TEST_SIMULATION_H
#define TEST_SIMULATION_H

#include "SimulationTime\SimulationTime.h"

#include <memory>

class SimulationParameter;
class SimulationResults;
class SimulationObserver;
class DataWriter;

class TestSimulation
{
public:
	//beachte Interface Segregation

	virtual std::shared_ptr< SimulationParameter> getSimulationParameter() = 0;
	virtual std::shared_ptr< DataWriter> getDataWriter() = 0;	
	
	virtual std::shared_ptr< SimulationResults> getSimulationResults() = 0;
	virtual bool getSimulationRun() = 0;
	virtual std::string getSimulationRunTimeOutput() = 0;
	virtual void registerSimulationObserver(std::shared_ptr< SimulationObserver> simObserver) = 0;
	
	
	virtual void makeSimulationHeadOutput() = 0;
	virtual void setStartTime() = 0;
	virtual void setEndTime() = 0;

private:

};
#endif