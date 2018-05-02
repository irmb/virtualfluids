#include "TestConditionFactoryImp.h"

#include "Utilities\TestParameter\TestParameter.h"
#include "Utilities\TestCondition\TestConditionImp.h"

std::shared_ptr<TestConditionFactory> TestConditionFactoryImp::getNewInstance(std::vector<std::shared_ptr<TestParameter>> testPara)
{
	return std::shared_ptr<TestConditionFactory>(new TestConditionFactoryImp(testPara));
}

TestConditionFactoryImp::TestConditionFactoryImp(std::vector<std::shared_ptr<TestParameter>> testPara):testPara(testPara)
{
}

std::vector<std::shared_ptr<TestCondition>> TestConditionFactoryImp::makeTestConditions()
{
	std::vector<std::shared_ptr<TestCondition>> testConditions;

	for (int i = 0; i < testPara.size(); i++) {
		std::shared_ptr<TestConditionImp> testCondit = TestConditionImp::getNewInstance();
		testCondit->initParameter(testPara.at(i)->getViscosity(), testPara.at(i)->getGridPath(), testPara.at(i)->getFilePath(), testPara.at(i)->getNumberOfGridLevels(), testPara.at(i)->getEndTime(), testPara.at(i)->getTimeStepLength());
		testCondit->initInitialConditions(testPara.at(i)->getInitialCondition());
		testCondit->initGridProvider();
		testCondit->initResults(testPara.at(i)->getLx(), testPara.at(i)->getLz(), testPara.at(i)->getTimeStepLength());
		testCondit->initDataWriter(testPara.at(i)->getYSliceForCalculation(), testPara.at(i)->getStartTimeCalculation(), testPara.at(i)->getEndTime(), testPara.at(i)->getTimeStepLength(), testPara.at(i)->getWriteFiles(), testPara.at(i)->getStartTimeDataWriter());
		testConditions.push_back(testCondit);
	}

	return testConditions;
}
