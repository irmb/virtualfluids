#include "TestConditionFactoryImp.h"

#include "Utilities\SimulationParameter\SimulationParameter.h"
#include "Utilities\TestCondition\TestConditionImp.h"

std::shared_ptr<TestConditionFactory> TestConditionFactoryImp::getNewInstance()
{
	return std::shared_ptr<TestConditionFactory>(new TestConditionFactoryImp());
}

TestConditionFactoryImp::TestConditionFactoryImp()
{

}

std::vector<std::shared_ptr<TestCondition>> TestConditionFactoryImp::makeTestConditions(std::vector< std::shared_ptr< SimulationParameter> > simPara)
{
	std::vector< std::shared_ptr<TestCondition> > testConditions;

	for (int i = 0; i < simPara.size(); i++) {
		std::shared_ptr< TestConditionImp> testCondit = TestConditionImp::getNewInstance();

		testCondit->initParameter(simPara.at(i)->getViscosity(), simPara.at(i)->getGridPath(), simPara.at(i)->getFilePath(), simPara.at(i)->getNumberOfGridLevels(), simPara.at(i)->getEndTime(), simPara.at(i)->getTimeStepLength(), simPara.at(i)->getDevices(), simPara.at(i)->getMaxVelocity());
		testCondit->initInitialConditions(simPara.at(i)->getInitialCondition());
		testCondit->initGridProvider();
		testCondit->initCalculator(simPara.at(i)->getCalculator());
		testCondit->initSimulationResults(simPara.at(i)->getLx(), simPara.at(i)->getLz(), simPara.at(i)->getTimeStepLength());
		testCondit->setTestResults(simPara.at(i)->getTestResults());
		testCondit->initDataWriter(simPara.at(i)->getYSliceForCalculation(), simPara.at(i)->getStartTimeCalculation(), simPara.at(i)->getEndTime(), simPara.at(i)->getTimeStepLength(), simPara.at(i)->getWriteFiles(), simPara.at(i)->getStartTimeDataWriter());
		testConditions.push_back(testCondit);
	}

	return testConditions;
}
