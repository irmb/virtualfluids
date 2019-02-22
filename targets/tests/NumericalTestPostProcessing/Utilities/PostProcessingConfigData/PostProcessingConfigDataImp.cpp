#include "PostProcessingConfigDataImp.h"

PostProcessingConfigDataImp::PostProcessingConfigDataImp()
{
}

std::shared_ptr<PostProcessingConfigDataImp> PostProcessingConfigDataImp::getNewInstance()
{
	return std::shared_ptr<PostProcessingConfigDataImp>(new PostProcessingConfigDataImp());
}

std::vector<BasicSimulation> PostProcessingConfigDataImp::getSimulations()
{
	return simulations;
}

std::vector<Assistant> PostProcessingConfigDataImp::getAssistants()
{
	return assistants;
}

std::vector<DataCombination> PostProcessingConfigDataImp::getDataCombinations()
{
	return dataCombinations;
}

void PostProcessingConfigDataImp::setSimulations(std::vector<BasicSimulation> simulations)
{
	this->simulations = simulations;
}

void PostProcessingConfigDataImp::setAssistants(std::vector<Assistant> assis)
{
	this->assistants = assis;
}

void PostProcessingConfigDataImp::setDataCombinations(std::vector<DataCombination> dataComb)
{
	this->dataCombinations = dataComb;
}
