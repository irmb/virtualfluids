#include "reader.h"

#include <fstream>
#include <iostream>

#include "utilities/input/Input.h"
#include "utilities/StringUtil/StringUtil.h"

#include "Utilities\EvaluationParameter\EvaluationParameter.h"
#include "Utilities\TestInformation\TestInformationImp.h"
#include "Utilities\TestParameter\TaylorGreenTestParameter\TaylorGreenTestParameter.h"
#include "Utilities\TestParameter\ShearWaveTestParameter\ShearWaveTestParameter.h"
#include "Utilities\TestResults\PhiAndNuTestResults.h"

void Reader::calcNumberOfEqualTests()
{
	for (int i = 0; i < tgv.size(); i++)
		if (tgv.at(i))
			numberOfTaylorGreenTests++;

	for (int i = 0; i < sw.size(); i++)
		if (sw.at(i))
			numberOfShearWaveTests++;
}


std::shared_ptr<Reader> Reader::getNewInstance(const std::string aFilePath)
{
	return std::shared_ptr<Reader>(new Reader(aFilePath));
}

Reader::Reader(const std::string aFilePath)
{
	std::ifstream stream;
	stream.open(aFilePath.c_str(), std::ios::in);
	std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

	viscosity = StringUtil::toDouble(input->getValue("Viscosity"));
	minOrderOfAccuracy = StringUtil::toDouble(input->getValue("MinOrderOfAccuracy"));

	amplitudeTGV = StringUtil::toDouble(input->getValue("Amplitude_TGV"));
	u0TGV = StringUtil::toDouble(input->getValue("u0_TGV"));
	v0SW = StringUtil::toDouble(input->getValue("v0_SW"));
	u0SW = StringUtil::toDouble(input->getValue("u0_SW"));
	
	numberOfTimeSteps = StringUtil::toInt(input->getValue("NumberOfTimeSteps"));
	basisTimeStepLength = StringUtil::toInt(input->getValue("BasisTimeStepLength"));
	startStepCalculation = StringUtil::toInt(input->getValue("StartStepCalculation"));

	grids.resize(5);
	grids.at(0) = input->getValue("GridPath32");
	grids.at(1) = input->getValue("GridPath64");
	grids.at(2) = input->getValue("GridPath128");
	grids.at(3) = input->getValue("GridPath256");
	grids.at(4) = input->getValue("GridPath512");

	ySliceForCalculation = StringUtil::toInt(input->getValue("ySliceForCalculation"));
	
	writeFiles = StringUtil::toBool(input->getValue("WriteFiles"));
	filePath = input->getValue("PathForFileWriting");
	startStepFileWriter = StringUtil::toInt(input->getValue("StartStepFileWriter"));
	logFilePath= input->getValue("PathLogFile");;

	tgv.resize(5);
	tgv.at(0) = StringUtil::toBool(input->getValue("TaylorGreenVortex32"));
	tgv.at(1) = StringUtil::toBool(input->getValue("TaylorGreenVortex64"));
	tgv.at(2) = StringUtil::toBool(input->getValue("TaylorGreenVortex128"));
	tgv.at(3) = StringUtil::toBool(input->getValue("TaylorGreenVortex256"));
	tgv.at(4) = StringUtil::toBool(input->getValue("TaylorGreenVortex512"));

	sw.resize(5);
	sw.at(0) = StringUtil::toBool(input->getValue("ShearWave32"));
	sw.at(1) = StringUtil::toBool(input->getValue("ShearWave64"));
	sw.at(2) = StringUtil::toBool(input->getValue("ShearWave128"));
	sw.at(3) = StringUtil::toBool(input->getValue("ShearWave256"));
	sw.at(4) = StringUtil::toBool(input->getValue("ShearWave512"));

	l.resize(5);
	l.at(0) = 32.0;
	l.at(1) = 64.0;
	l.at(2) = 128.0;
	l.at(3) = 256.0;
	l.at(4) = 512.0;

	stream.close();

	numberOfTaylorGreenTests = 0;
	numberOfShearWaveTests = 0;
	calcNumberOfEqualTests();
}

std::vector<std::shared_ptr<EvaluationParameter>> Reader::makeEvaluationParameter()
{
	std::vector<std::shared_ptr<EvaluationParameter>> evaPara;

	for (int i = 0; i < tgv.size(); i++) {
		if (tgv.at(i)) {
			evaPara.push_back(EvaluationParameter::getNewInstance("TaylorGreenVortex", numberOfTaylorGreenTests, l.at(i),"vX", logFilePath, minOrderOfAccuracy, writeFiles, viscosity));
		}
	}
	for (int i = 0; i < sw.size(); i++) {
		if (sw.at(i)) {
			evaPara.push_back(EvaluationParameter::getNewInstance("ShearWave", numberOfShearWaveTests, l.at(i),"vZ", logFilePath, minOrderOfAccuracy, writeFiles, viscosity));
		}
	}
	return evaPara;
}

std::shared_ptr<TestInformation> Reader::makeTestInformation()
{
	bool tgvTest = false;
	bool swTest = false;
	for (int i = 0; i < tgv.size(); i++) {
		if (tgv.at(i))
			tgvTest = true;
	}
	for (int i = 0; i < sw.size(); i++) {
		if (sw.at(i))
			swTest = true;
	}
	std::shared_ptr<TestInformation> testInfo = TestInformationImp::getNewInstance(numberOfTimeSteps, basisTimeStepLength, startStepCalculation, viscosity, tgvTest, u0TGV, amplitudeTGV, swTest, u0SW, v0SW);
	return testInfo;
}

std::vector<std::shared_ptr<TestParameter>> Reader::makeTestParameter()
{
	std::vector<std::shared_ptr<TestParameter>> testParameter;
	std::shared_ptr<PhiAndNuTestResults> tgvTestResults = PhiAndNuTestResults::getNewInstance("TaylorGreenVortex");
	std::shared_ptr<PhiAndNuTestResults> swTestResults = PhiAndNuTestResults::getNewInstance("ShearWave");

	for (int i = 0; i < tgv.size(); i++) {
		if (tgv.at(i)) {
			testParameter.push_back(TaylorGreenTestParameter::getNewInstance(u0TGV, amplitudeTGV, viscosity, l.at(i), numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, grids.at(i), writeFiles, startStepFileWriter, filePath, tgvTestResults));
		}
	}
	for (int i = 0; i < sw.size(); i++) {
		if (sw.at(i)) {
			testParameter.push_back(ShearWaveTestParameter::getNewInstance(u0SW, v0SW, viscosity, l.at(i), numberOfTimeSteps, basisTimeStepLength, startStepCalculation, ySliceForCalculation, grids.at(i), writeFiles, startStepFileWriter, filePath, swTestResults));
		}
	}

	return testParameter;
}
