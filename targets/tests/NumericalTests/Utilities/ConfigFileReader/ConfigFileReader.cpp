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

	configData->l0 = 32.0;
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
	configData->viscosity = StringUtil::toDoubleVector(input->getValue("Viscosity"));
	configData->minOrderOfAccuracy = StringUtil::toDouble(input->getValue("MinOrderOfAccuracy"));
	configData->dataToCalcPhiAndNuTest = StringUtil::toString(input->getValue("DataToCalc_PhiAndNu"));
	configData->startTimeStepCalculationPhiNu = StringUtil::toInt(input->getValue("StartTimeStepCalculation_PhiNu"));
	configData->endTimeStepCalculationPhiNu = StringUtil::toInt(input->getValue("EndTimeStepCalculation_PhiNu"));
	configData->maxL2NormDiff = StringUtil::toDouble(input->getValue("MaxL2NormDiff"));
	configData->dataToCalcL2Test = StringUtil::toString(input->getValue("DataToCalc_L2"));
	configData->basicTimeStepL2Norm = StringUtil::toInt(input->getValue("BasicTimeStep_L2"));
	configData->divergentTimeStepL2Norm = StringUtil::toInt(input->getValue("DivergentTimeStep_L2"));
	configData->amplitudeTGV = StringUtil::toDoubleVector(input->getValue("Amplitude_TGV"));
	configData->u0TGV = StringUtil::toDoubleVector(input->getValue("u0_TGV"));
	configData->nuAndPhiTestTGV = StringUtil::toBool(input->getValue("PhiAndNuTest_TGV"));
	configData->l2NormTestTGV = StringUtil::toBool(input->getValue("L2NormTest_TGV"));
	configData->v0SW = StringUtil::toDoubleVector(input->getValue("v0_SW"));
	configData->u0SW = StringUtil::toDoubleVector(input->getValue("u0_SW"));
	configData->nuAndPhiTestSW = StringUtil::toBool(input->getValue("PhiAndNuTest_SW"));
	configData->l2NormTestSW = StringUtil::toBool(input->getValue("L2NormTest_SW"));
	configData->numberOfTimeSteps = StringUtil::toInt(input->getValue("NumberOfTimeSteps"));
	configData->basisTimeStepLength = StringUtil::toInt(input->getValue("BasisTimeStepLength"));
	configData->grids.resize(5);
	configData->grids.at(0) = input->getValue("GridPath32");
	configData->grids.at(1) = input->getValue("GridPath64");
	configData->grids.at(2) = input->getValue("GridPath128");
	configData->grids.at(3) = input->getValue("GridPath256");
	configData->grids.at(4) = input->getValue("GridPath512");
	configData->numberOfGridLevels = StringUtil::toInt(input->getValue("NumberOfGridLevels"));
	configData->maxLevel = configData->numberOfGridLevels - 1;
	configData->ySliceForCalculation = StringUtil::toInt(input->getValue("ySliceForCalculation"));
	configData->writeFiles = StringUtil::toBool(input->getValue("WriteFiles"));
	configData->filePath = input->getValue("PathForFileWriting");
	configData->startStepFileWriter = StringUtil::toInt(input->getValue("StartStepFileWriter"));
	configData->logFilePath = input->getValue("PathLogFile");;
	configData->tgv.resize(5);
	configData->tgv.at(0) = StringUtil::toBool(input->getValue("TaylorGreenVortex32"));
	configData->tgv.at(1) = StringUtil::toBool(input->getValue("TaylorGreenVortex64"));
	configData->tgv.at(2) = StringUtil::toBool(input->getValue("TaylorGreenVortex128"));
	configData->tgv.at(3) = StringUtil::toBool(input->getValue("TaylorGreenVortex256"));
	configData->tgv.at(4) = StringUtil::toBool(input->getValue("TaylorGreenVortex512"));
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
	if (configData->u0TGV.size() != configData->amplitudeTGV.size()) {
		std::cout << "Length u0_TGV is unequal to Lenght Amplitude_TGV!" << std::endl << std::flush;
		exit(1);
	}
		
	if (configData->u0SW.size() != configData->v0SW.size()) {
		std::cout << "Length u0_SW is unequal to Lenght v0_SW!" << std::endl << std::flush;
		exit(1);
	}	
}