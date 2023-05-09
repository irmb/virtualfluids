#include "PostProcessingConfigFileReaderImp.h"

#include <basics/config/ConfigurationFile.h>
#include "StringUtilities/StringUtil.h"

#include "Utilities/PostProcessingConfigData/PostProcessingConfigDataImp.h"

#include <fstream>


std::shared_ptr<PostProcessingConfigFileReader> PostProcessingConfigFileReaderImp::getNewInstance()
{
	return std::shared_ptr<PostProcessingConfigFileReader>(new PostProcessingConfigFileReaderImp());
}

std::shared_ptr<PostProcessingConfigData> PostProcessingConfigFileReaderImp::readConfigFile(std::string filePath)
{
	auto input = std::make_shared<vf::basics::ConfigurationFile>();
	input->load(filePath);

	std::vector<BasicSimulation> simulation;
	std::vector<Assistant> assistants;
	std::vector<DataCombination> combination;

	if(StringUtil::toBool(input->getValue<std::string>("ShearWave")))
		simulation.push_back(ShearWave);

	if (StringUtil::toBool(input->getValue<std::string>("TaylorGreenVortexUx")))
		simulation.push_back(TaylorGreenVortexUx);

	if (StringUtil::toBool(input->getValue<std::string>("TaylorGreenVortexUz")))
		simulation.push_back(TaylorGreenVortexUz);

	if (StringUtil::toBool(input->getValue<std::string>("Phi")))
		assistants.push_back(Phi);

	if (StringUtil::toBool(input->getValue<std::string>("Ny")))
		assistants.push_back(Ny);

	if (StringUtil::toBool(input->getValue<std::string>("L2Norm")))
		assistants.push_back(L2Norm);

	if (StringUtil::toBool(input->getValue<std::string>("L2Norm_BetweenKernels")))
		assistants.push_back(L2NormBetweenKernels);

	if (StringUtil::toBool(input->getValue<std::string>("TimeOutput")))
		assistants.push_back(Time);


	if (StringUtil::toBool(input->getValue<std::string>("EqualSimulationsForDifferentKernels")))
		combination.push_back(EqualSimulationsForDifferentKernels);

	if (StringUtil::toBool(input->getValue<std::string>("EqualKernelSimulationsForDifferentViscosities")))
		combination.push_back(EqualKernelSimulationsForDifferentViscosities);

	std::shared_ptr<PostProcessingConfigDataImp> data = PostProcessingConfigDataImp::getNewInstance();

	data->setAssistants(assistants);
	data->setSimulations(simulation);
	data->setDataCombinations(combination);

	data->setLogFilesPath(input->getValue<std::string>("LogFilesPath"));
	data->setMathematicaFilePath(input->getValue<std::string>("MathematicaFilePath"));
	
	return data;
}

PostProcessingConfigFileReaderImp::PostProcessingConfigFileReaderImp()
{
}
