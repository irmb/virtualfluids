#ifndef NUMERICAL_TEST_SIMULATION_H
#define NUMERICAL_TEST_SIMULATION_H

#include <memory>

class SimulationObserver; 

enum SimulationStatus { initialized , executed, crashed};

class NumericalTestSimulation
{
public:
	virtual ~NumericalTestSimulation() = default;
	virtual void run() = 0;
	virtual SimulationStatus getSimulationStatus() = 0;
	virtual void registerSimulationObserver(std::shared_ptr<SimulationObserver> simObserver) = 0;
};
#endif