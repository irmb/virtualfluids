#include "ConfigFileReader.h"

#include "utilities/input/Input.h"
#include "utilities/StringUtil/StringUtil.h"

#include <fstream>
#include <iostream>

std::shared_ptr<ConfigFileReader> ConfigFileReader::getNewInstance(const std::string aFilePath)
{
	return std::shared_ptr<ConfigFileReader>(new ConfigFileReader(aFilePath));
}

ConfigFileReader::ConfigFileReader(const std::string aFilePath)
{
	configData = std::shared_ptr< ConfigDataStruct>(new ConfigDataStruct);

	configData->lx.resize(5);
	configData->lx.at(0) = 32.0;
	configData->lx.at(1) = 64.0;
	configData->lx.at(2) = 128.0;
	configData->lx.at(3) = 256.0;
	configData->lx.at(4) = 512.0;

	configData->lz.resize(5);
	configData->lz.at(0) = configData->lx.at(0) * 3.0 / 2.0;
	configData->lz.at(1) = configData->lx.at(1) * 3.0 / 2.0;
	configData->lz.at(2) = configData->lx.at(2) * 3.0 / 2.0;
	configData->lz.at(3) = configData->lx.at(3) * 3.0 / 2.0;
	configData->lz.at(4) = configData->lx.at(4) * 3.0 / 2.0;

	configData->rho0 = 1.0;

	readConfigFile(aFilePath);
}

void ConfigFileReader::readConfigFile(const std::string aFilePath)
{
	std::ifstream stream;
	stream.open(aFilePath.c_str(), std::ios::in);
	if (stream.fail())
		throw "can not open config file!\n";

	std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

	configData->devices = StringUtil::toIntVector(input->getValue("Devices"));
	configData->kernelsToTest = StringUtil::toStringVector(input->getValue("KernelsToTest"));
	configData->numberOfTimeSteps = StringUtil::toInt(input->getValue("NumberOfTimeSteps"));
	configData->viscosity = StringUtil::toDoubleVector(input->getValue("Viscosity"));
	configData->minOrderOfAccuracy = StringUtil::toDouble(input->getValue("MinOrderOfAccuracy"));
	configData->dataToCalcPhiAndNuTest = StringUtil::toStringVector(input->getValue("DataToCalc_PhiAndNu"));
	configData->startTimeStepCalculationPhiNu = StringUtil::toInt(input->getValue("StartTimeStepCalculation_PhiNu"));
	configData->endTimeStepCalculationPhiNu = StringUtil::toInt(input->getValue("EndTimeStepCalculation_PhiNu"));
	configData->nuAndPhiTest = StringUtil::toBool(input->getValue("PhiAndNuTest"));
	configData->l2NormTest = StringUtil::toBool(input->getValue("L2NormTest"));
	configData->maxL2NormDiff = StringUtil::toDouble(input->getValue("MaxL2NormDiff"));
	configData->dataToCalcL2Test = StringUtil::toStringVector(input->getValue("DataToCalc_L2"));
	configData->basicTimeStepL2Norm = StringUtil::toInt(input->getValue("BasicTimeStep_L2"));
	configData->divergentTimeStepL2Norm = StringUtil::toInt(input->getValue("DivergentTimeStep_L2"));
	configData->basisTimeStepLengthTGVux = StringUtil::toIntVector(input->getValue("BasisTimeStepLength_TGV_Ux"));
	configData->amplitudeTGVux = StringUtil::toDoubleVector(input->getValue("Amplitude_TGV_Ux"));
	configData->u0TGVux = StringUtil::toDoubleVector(input->getValue("ux_TGV_Ux"));
	configData->l0TGVux = StringUtil::toInt(input->getValue("l0_TGV_Ux"));
	configData->basisTimeStepLengthTGVuz = StringUtil::toIntVector(input->getValue("BasisTimeStepLength_TGV_Uz"));
	configData->amplitudeTGVuz = StringUtil::toDoubleVector(input->getValue("Amplitude_TGV_Uz"));
	configData->v0TGVuz = StringUtil::toDoubleVector(input->getValue("uz_TGV_Uz"));
	configData->l0TGVuz = StringUtil::toInt(input->getValue("l0_TGV_Uz"));

	configData->basisTimeStepLengthSW = StringUtil::toIntVector(input->getValue("BasisTimeStepLength_SW"));
	configData->v0SW = StringUtil::toDoubleVector(input->getValue("v0_SW"));
	configData->u0SW = StringUtil::toDoubleVector(input->getValue("u0_SW"));
	configData->l0SW = StringUtil::toInt(input->getValue("l0_SW"));

	configData->l2NormBetweenKernelTest = StringUtil::toBool(input->getValue("L2NormBetweenKernelsTest"));
	configData->basicKernelL2NormTest = StringUtil::toString(input->getValue("BasicKernel_L2NormBetweenKernels"));
	configData->timeStepsL2NormBetweenKernel = StringUtil::toIntVector(input->getValue("Timesteps_L2NormBetweenKernels"));
	configData->dataToCalcL2NormBetweenKernel = StringUtil::toStringVector(input->getValue("DataToCalc_L2NormBetweenKernels"));
	
	configData->grids.resize(5);
	configData->grids.at(0) = input->getValue("GridPath32");
	configData->grids.at(1) = input->getValue("GridPath64");
	configData->grids.at(2) = input->getValue("GridPath128");
	configData->grids.at(3) = input->getValue("GridPath256");
	configData->grids.at(4) = input->getValue("GridPath512");
	configData->numberOfGridLevels = StringUtil::toInt(input->getValue("NumberOfGridLevels"));
	configData->maxLevel = configData->numberOfGridLevels - 1;
	configData->ySliceForCalculation = StringUtil::toInt(input->getValue("ySliceForCalculation"));
	configData->writeFiles = StringUtil::toBool(input->getValue("WriteVTKFiles"));
	configData->writeAnalyticalToVTK = StringUtil::toBool(input->getValue("WriteAnalyResultsToVTK"));
	configData->filePath = input->getValue("PathForFileWriting");
	configData->startStepFileWriter = StringUtil::toInt(input->getValue("StartStepFileWriter"));
	configData->logFilePath = input->getValue("PathLogFile");
	configData->tgvUx.resize(5);
	configData->tgvUx.at(0) = StringUtil::toBool(input->getValue("TaylorGreenVortexUx32"));
	configData->tgvUx.at(1) = StringUtil::toBool(input->getValue("TaylorGreenVortexUx64"));
	configData->tgvUx.at(2) = StringUtil::toBool(input->getValue("TaylorGreenVortexUx128"));
	configData->tgvUx.at(3) = StringUtil::toBool(input->getValue("TaylorGreenVortexUx256"));
	configData->tgvUx.at(4) = StringUtil::toBool(input->getValue("TaylorGreenVortexUx512"));
	configData->tgvUz.resize(5);
	configData->tgvUz.at(0) = StringUtil::toBool(input->getValue("TaylorGreenVortexUz32"));
	configData->tgvUz.at(1) = StringUtil::toBool(input->getValue("TaylorGreenVortexUz64"));
	configData->tgvUz.at(2) = StringUtil::toBool(input->getValue("TaylorGreenVortexUz128"));
	configData->tgvUz.at(3) = StringUtil::toBool(input->getValue("TaylorGreenVortexUz256"));
	configData->tgvUz.at(4) = StringUtil::toBool(input->getValue("TaylorGreenVortexUz512"));
	configData->sw.resize(5);
	configData->sw.at(0) = StringUtil::toBool(input->getValue("ShearWave32"));
	configData->sw.at(1) = StringUtil::toBool(input->getValue("ShearWave64"));
	configData->sw.at(2) = StringUtil::toBool(input->getValue("ShearWave128"));
	configData->sw.at(3) = StringUtil::toBool(input->getValue("ShearWave256"));
	configData->sw.at(4) = StringUtil::toBool(input->getValue("ShearWave512"));

	stream.close();

	checkConfigFileData();
}

std::shared_ptr<ConfigDataStruct> ConfigFileReader::getConfigData()
{
	return configData;
}

void ConfigFileReader::checkConfigFileData()
{
	if (configData->u0TGVux.size() != configData->amplitudeTGVux.size() || configData->u0TGVux.size() != configData->basisTimeStepLengthTGVux.size()) {
		std::cout << "Length u0_TGV_U0 is unequal to Lenght Amplitude_TGV_U0 or BasisTimeStepLength_TGV_U0!" << std::endl << std::flush;
		exit(1);
	}

	if (configData->v0TGVuz.size() != configData->amplitudeTGVuz.size() || configData->v0TGVuz.size() != configData->basisTimeStepLengthTGVuz.size()) {
		std::cout << "Length v0_TGV_V0 is unequal to Lenght Amplitude_TGV_V0 or BasisTimeStepLength_TGV_V0!" << std::endl << std::flush;
		exit(1);
	}
		
	if (configData->u0SW.size() != configData->v0SW.size() || configData->u0SW.size() != configData->basisTimeStepLengthSW.size()) {
		std::cout << "Length u0_SW is unequal to Lenght v0_SW!" << std::endl << std::flush;
		exit(1);
	}	
}