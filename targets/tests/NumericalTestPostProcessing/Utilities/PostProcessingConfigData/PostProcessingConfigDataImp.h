#ifndef POST_PROCESSING_CONFIG_DATA_IMP_H
#define POST_PROCESSING_CONFIG_DATA_IMP_H

#include "PostProcessingConfigData.h"

class PostProcessingConfigDataImp : public PostProcessingConfigData
{
public:
	static std::shared_ptr<PostProcessingConfigDataImp> getNewInstance();

	std::vector<BasicSimulation> getSimulations();
	std::vector<Assistant> getAssistants();
	std::vector<DataCombination> getDataCombinations();

	void setSimulations(std::vector<BasicSimulation> simulations);
	void setAssistants(std::vector<Assistant> assis);
	void setDataCombinations(std::vector<DataCombination> dataComb);

private:
	PostProcessingConfigDataImp();

	std::vector<BasicSimulation> simulations;
	std::vector<Assistant> assistants;
	std::vector<DataCombination> dataCombinations;

};
#endif