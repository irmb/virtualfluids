#ifndef TEST_CONDITION_IMP_H
#define TEST_CONDITION_IMP_H

#include "TestCondition.h"

#include <string>
#include <memory>

class InitialCondition;
class FileWriter;
class ToVectorWriter;
class TestParameter;

class TestConditionImp: public TestCondition
{
public:
	static std::shared_ptr<TestConditionImp> getNewInstance();
	void initParameter(real viscosity, std::string gridPath, std::string filePath, int numberOfGridLevels, unsigned int endTime, unsigned int timeStepLength);
	void initInitialConditions(std::shared_ptr<InitialCondition> initialCondition);
	void initGridProvider();
	void initResults(unsigned int lx, unsigned int lz, unsigned int timeStepLength);
	void initDataWriter(unsigned int ySliceForCalculation, unsigned int startTimeCalculation, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, unsigned int startTimeDataWriter);

	std::shared_ptr<Parameter> getParameter();
	std::shared_ptr<GridProvider> getGrid();
	std::shared_ptr<DataWriter> getDataWriter();
	std::shared_ptr<Results> getSimulationResults();

protected:
	TestConditionImp() {};
		
private:
	std::shared_ptr<TestParameter> testPara;
	std::shared_ptr<Parameter> para;
	std::shared_ptr<InitialCondition> initialCondition;
	std::shared_ptr<GridProvider> grid;
	std::shared_ptr<Results> simResults;
	std::shared_ptr<ToVectorWriter> writeToVector;
	std::shared_ptr<FileWriter> fileWriter;
	
};
#endif