#ifndef POST_PROCESSING_STRATEGY_H
#define POST_PROCESSING_STRATEGY_H

#include <vector>
#include <memory>
#include <string>

class SimulationResults;

class PostProcessingStrategy
{
public:
	virtual void evaluate() = 0;
};
#endif