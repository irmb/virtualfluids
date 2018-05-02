#ifndef READER_H
#define READER_H

#include "LBM\LB.h"

#include <memory>
#include <vector>
#include <string>

class EvaluationParameter;
class TestInformation;
class TestParameter;

class Reader
{
public:
	static std::shared_ptr < Reader > getNewInstance(const std::string aFilePath);
	std::vector < std::shared_ptr < EvaluationParameter > > makeEvaluationParameter();
	std::shared_ptr < TestInformation > makeTestInformation();
	std::vector < std::shared_ptr < TestParameter > > makeTestParameter();

protected:
	Reader() {};
	Reader(const std::string aFilePath);
	
private:
	void calcNumberOfEqualTests();

	real viscosity;
	double minOrderOfAccuracy;

	real u0SW, v0SW;
	real amplitudeTGV, u0TGV;

	unsigned int numberOfTimeSteps, basisTimeStepLength, startStepCalculation;
	unsigned int ySliceForCalculation;
	std::vector<real> l;

	std::vector<std::string> grids;

	bool writeFiles;
	std::string filePath;
	unsigned int startStepFileWriter;
	std::string logFilePath;
	

	std::vector<bool> tgv;
	std::vector<bool> sw;

	int numberOfTaylorGreenTests, numberOfShearWaveTests;
};
#endif