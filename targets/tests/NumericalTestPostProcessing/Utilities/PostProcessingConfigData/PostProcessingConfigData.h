#ifndef POST_PROCESSING_CONFIG_DATA_H
#define POST_PROCESSING_CONFIG_DATA_H

#include "Simulation/BasicSimulation.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistant.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistant.h"

#include <string>
#include <vector>

class PostProcessingConfigData
{
public:
	virtual std::vector<BasicSimulation> getSimulations() = 0;
	virtual std::vector<Assistant> getAssistants() = 0;
	virtual std::vector<DataCombination> getDataCombinations() = 0;

	virtual std::string getMathematicaFilePath() = 0;
	virtual std::string getLogFilesPath() = 0;
};
#endif