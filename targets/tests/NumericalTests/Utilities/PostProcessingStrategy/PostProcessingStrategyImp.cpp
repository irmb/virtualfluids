#include "PostProcessingStrategyImp.h"

#include "Utilities\Results\SimulationResults\SimulationResults.h"

bool PostProcessingStrategyImp::checkEqualSimulationResults(std::shared_ptr<SimulationResults> simResultsToCheck)
{ 
	if (simResultsToCheck.get() == simResult.get())
		return true;

	return false;
}

int PostProcessingStrategyImp::getNumberOfXNodes()
{
	return simResult->getNumberOfXNodes();
}

int PostProcessingStrategyImp::getNumberOfYNodes()
{
	return simResult->getNumberOfYNodes();
}

int PostProcessingStrategyImp::getNumberOfZNodes()
{
	return simResult->getNumberOfZNodes();
}

PostProcessingStrategyImp::PostProcessingStrategyImp(std::shared_ptr<SimulationResults> simResult) : simResult(simResult)
{
}

int PostProcessingStrategyImp::calcTimeStepInResults(unsigned int timeStep)
{
	for (int i = 0; i < simResult->getTimeSteps().size(); i++) {
		if (timeStep == simResult->getTimeSteps().at(i))
			return simResult->getTimeSteps().at(i);
	}
}
