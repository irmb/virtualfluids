#ifndef TEST_CONDITION_IMP_H
#define TEST_CONDITION_IMP_H

#include "TestCondition.h"

#include <string>
#include <memory>

class InitialCondition;
class FileWriter;
class ToVectorWriter;
class TestParameter;
class Calculator;
class TestResults;
class SimulationResults;

class TestConditionImp: public TestCondition
{
public:
	static std::shared_ptr<TestConditionImp> getNewInstance();
	void initParameter(real viscosity, std::string gridPath, std::string filePath, int numberOfGridLevels, unsigned int endTime, unsigned int timeStepLength, std::vector<int> devices, real velocity);
	void initInitialConditions(std::shared_ptr<InitialCondition> initialCondition);
	void initGridProvider();
	void initCalculator(std::shared_ptr<Calculator> calc);
	void initDataWriter(unsigned int ySliceForCalculation, unsigned int startTimeCalculation, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, unsigned int startTimeDataWriter);
	void initSimulationResults(unsigned int lx, unsigned int lz, unsigned int timeStepLength);

	void setTestResults(std::shared_ptr<TestResults> testResults);

	std::shared_ptr<Parameter> getParameter();
	std::shared_ptr<GridProvider> getGrid();
	std::shared_ptr<DataWriter> getDataWriter();
	std::shared_ptr<TestResults> getTestResults();
	std::shared_ptr<Calculator> getCalculator();

protected:
	TestConditionImp() {};
		
private:
	std::shared_ptr<TestParameter> testPara;
	std::shared_ptr<Parameter> para;
	std::shared_ptr<InitialCondition> initialCondition;
	std::shared_ptr<GridProvider> grid;
	std::shared_ptr<SimulationResults> simResults;
	std::shared_ptr<ToVectorWriter> writeToVector;
	std::shared_ptr<FileWriter> fileWriter;
	std::shared_ptr<Calculator> calculator;
	std::shared_ptr<TestResults> testResults;	
};
#endif