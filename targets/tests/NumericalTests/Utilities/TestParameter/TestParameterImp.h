#ifndef TEST_PARAMETER_IMP_H
#define TEST_PARAMETER_IMP_H

#include "TestParameter.h"

#include "LBM\LB.h"

class SimulationResults;

class TestParameterImp: public TestParameter
{
public:
	double getViscosity();
	std::string getGridPath();
	std::string getFilePath();
	unsigned int getNumberOfGridLevels();
	unsigned int getEndTime();
	unsigned int getTimeStepLength();
	unsigned int getLx();
	unsigned int getLz();
	unsigned int getYSliceForCalculation();
	unsigned int getStartTimeCalculation();
	bool getWriteFiles();
	unsigned int getStartTimeDataWriter();
	std::vector<int> getDevices();
	std::shared_ptr<InitialCondition> getInitialCondition();
	std::shared_ptr<Calculator> getCalculator();
	std::shared_ptr<TestResults> getTestResults();

protected:
	TestParameterImp() {};
	TestParameterImp(real viscosity, unsigned int lx, unsigned int lz, unsigned int l0,
		unsigned int numberOfTimeSteps, unsigned int basisTimeStepLength,
		unsigned int startStepCalculation, unsigned int ySliceForCalculation,
		std::string gridPath, unsigned int maxLevel, unsigned int numberOfGridLevels,
		bool writeFiles, unsigned int startStepFileWriter,
		std::shared_ptr<TestResults> testResults,
		std::vector<int> devices);

	std::string filePath;
	real viscosity;
	unsigned int lx;
	unsigned int numberOfTimeSteps, basisTimeStepLength;
	unsigned int startStepCalculation, startStepFileWriter, ySliceForCalculation;
	std::string gridPath;
	bool writeFiles;
	std::vector<int> devices;

	unsigned int maxLevel, numberOfGridLevels;
	unsigned int l0, lz;
	unsigned int timeStepLength;
	unsigned int startTimeCalculation, startTimeDataWriter;
	unsigned int endTime;

	std::shared_ptr<InitialCondition> initialCondition;
	std::shared_ptr<Calculator> calculator;
	std::shared_ptr<SimulationResults> simResults;
	std::shared_ptr<TestResults> testResults;
};

#endif
