#include "PostProcessingConfigFileReaderImp.h"

#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"

#include "Utilities/PostProcessingConfigData/PostProcessingConfigDataImp.h"

#include <fstream>


std::shared_ptr<PostProcessingConfigFileReader> PostProcessingConfigFileReaderImp::getNewInstance()
{
	return std::shared_ptr<PostProcessingConfigFileReader>(new PostProcessingConfigFileReaderImp());
}

std::shared_ptr<PostProcessingConfigData> PostProcessingConfigFileReaderImp::readConfigFile(std::string filePath)
{
	std::ifstream stream;
	stream.open(filePath.c_str(), std::ios::in);
	if (stream.fail()) {
		throw "can not open config file!\n";
		exit(1);
	}
	std::shared_ptr<input::Input> input = input::Input::makeInput(stream, "config");


	std::vector<BasicSimulation> simulation;
	std::vector<Assistant> assistants;
	std::vector<DataCombination> combination;

	if(StringUtil::toBool(input->getValue("ShearWave")))
		simulation.push_back(ShearWave);

	if (StringUtil::toBool(input->getValue("TaylorGreenVortexUx")))
		simulation.push_back(TaylorGreenVortexUx);

	if (StringUtil::toBool(input->getValue("TaylorGreenVortexUz")))
		simulation.push_back(TaylorGreenVortexUz);

	if (StringUtil::toBool(input->getValue("Phi")))
		assistants.push_back(Phi);

	if (StringUtil::toBool(input->getValue("Ny")))
		assistants.push_back(Ny);

	if (StringUtil::toBool(input->getValue("L2Norm")))
		assistants.push_back(L2Norm);

	if (StringUtil::toBool(input->getValue("L2Norm_BetweenKernels")))
		assistants.push_back(L2NormBetweenKernels);

	if (StringUtil::toBool(input->getValue("TimeOutput")))
		assistants.push_back(Time);


	if (StringUtil::toBool(input->getValue("EqualSimulationsForDifferentKernels")))
		combination.push_back(EqualSimulationsForDifferentKernels);

	if (StringUtil::toBool(input->getValue("EqualKernelSimulationsForDifferentViscosities")))
		combination.push_back(EqualKernelSimulationsForDifferentViscosities);

	std::shared_ptr<PostProcessingConfigDataImp> data = PostProcessingConfigDataImp::getNewInstance();

	data->setAssistants(assistants);
	data->setSimulations(simulation);
	data->setDataCombinations(combination);

	data->setLogFilesPath(input->getValue("LogFilesPath"));
	data->setMathematicaFilePath(input->getValue("MathematicaFilePath"));
	
	return data;
}

PostProcessingConfigFileReaderImp::PostProcessingConfigFileReaderImp()
{
}
