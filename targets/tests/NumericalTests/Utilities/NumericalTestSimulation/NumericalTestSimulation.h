#ifndef NUMERICAL_TEST_SIMULATION_H
#define NUMERICAL_TEST_SIMULATION_H

class SimulationObserver; 

class NumericalTestSimulation
{
public:
	virtual bool getSimulationRun() = 0;
	virtual void registerSimulationObserver(std::shared_ptr<SimulationObserver> simObserver) = 0;
};
#endif