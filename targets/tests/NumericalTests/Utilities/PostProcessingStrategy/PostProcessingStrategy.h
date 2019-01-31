#ifndef POST_PROCESSING_STRATEGY_H
#define POST_PROCESSING_STRATEGY_H

#include <vector>
#include <memory>

class SimulationResults;

class PostProcessingStrategy
{
public:
	virtual void evaluate() = 0;
	virtual bool checkEqualSimulationResults(std::shared_ptr<SimulationResults> simResultsToCheck) = 0;

};
#endif