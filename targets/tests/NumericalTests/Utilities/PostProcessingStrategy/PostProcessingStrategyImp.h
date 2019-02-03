#ifndef POST_PROCESSING_STRATEGY_IMP_H
#define POST_PROCESSING_STRATEGY_IMP_H

#include "PostProcessingStrategy.h"

#include <memory>

class SimulationResults;

class PostProcessingStrategyImp : public PostProcessingStrategy
{
public:
	virtual void evaluate() = 0;

	int getNumberOfXNodes();
	int getNumberOfYNodes();
	int getNumberOfZNodes();

protected:
	PostProcessingStrategyImp(std::shared_ptr<SimulationResults> simResult);
	int calcTimeStepInResults(unsigned int timeStep);

	std::shared_ptr<SimulationResults> simResult;

private:
	PostProcessingStrategyImp() {};
};
#endif