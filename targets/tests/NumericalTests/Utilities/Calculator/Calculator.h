#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <memory>

class SimulationResults;
class TestResults;

class Calculator {
public:
	virtual void calcAndCopyToTestResults() = 0;
	virtual void setSimulationResults(std::shared_ptr<SimulationResults> simResults) = 0;
};
#endif