#ifndef CONFIG_FILE_READER_H
#define CONFIG_FILE_READER_H

#include "VirtualFluids_GPU/Kernel//Utilities/KernelType.h"

#include "Core/Input/Input.h"
#include "Utilities/Structs/ConfigDataStruct.h"

#include <memory>
#include <string>

class KernelMapper;

class ConfigFileReader
{
public:
	static std::shared_ptr<ConfigFileReader> getNewInstance(const std::string aFilePath);
	std::shared_ptr<ConfigDataStruct> getConfigData();
	void readConfigFile();
	
private:
	ConfigFileReader() {};
	ConfigFileReader(const std::string aFilePath);
	
	std::ifstream openConfigFile(const std::string aFilePath);
	bool checkConfigFile(std::shared_ptr<input::Input> input);
	std::vector<KernelType> readKernelList(std::shared_ptr<input::Input> input);

	int calcNumberOfSimulations(std::shared_ptr<input::Input> input);
	int calcNumberOfSimulationGroup(std::shared_ptr<input::Input> input, std::string simName);
	unsigned int calcStartStepForToVectorWriter(std::shared_ptr<input::Input> input);

	std::vector<std::shared_ptr<TaylorGreenVortexUxParameterStruct> > makeTaylorGreenVortexUxParameter(std::shared_ptr<input::Input> input, std::shared_ptr<BasicSimulationParameterStruct> basicSimParameter);
	std::vector<std::shared_ptr<TaylorGreenVortexUzParameterStruct> > makeTaylorGreenVortexUzParameter(std::shared_ptr<input::Input> input, std::shared_ptr<BasicSimulationParameterStruct> basicSimParameter);
	std::vector<std::shared_ptr<ShearWaveParameterStruct> > makeShearWaveParameter(std::shared_ptr<input::Input> input, std::shared_ptr<BasicSimulationParameterStruct> basicSimParameter);

	std::shared_ptr<NyTestParameterStruct> makeNyTestParameter(std::shared_ptr<input::Input> input);
	std::shared_ptr<PhiTestParameterStruct> makePhiTestParameter(std::shared_ptr<input::Input> input);
	std::shared_ptr<L2NormTestParameterStruct> makeL2NormTestParameter(std::shared_ptr<input::Input> input);
	std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> makeL2NormTestBetweenKernelsParameter(std::shared_ptr<input::Input> input);

	std::shared_ptr<BasicSimulationParameterStruct> makeBasicSimulationParameter(std::shared_ptr<input::Input> input);
	std::vector<std::shared_ptr<GridInformationStruct> > makeGridInformation(std::shared_ptr<input::Input> input, std::string simName);

	std::shared_ptr<VectorWriterInformationStruct> makeVectorWriterInformationStruct(std::shared_ptr<input::Input> input);
	std::shared_ptr<LogFileParameterStruct> makeLogFilePara(std::shared_ptr<input::Input> input);
	

	const std::string myFilePath;
	std::shared_ptr<ConfigDataStruct> configData;
	std::shared_ptr<KernelMapper> myKernelMapper;
};
#endif