#include "L2NormTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\Calculator\L2NormCalculator\L2NormCalculator.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"

std::shared_ptr<L2NormTest> L2NormTest::getNewInstance(std::shared_ptr< AnalyticalResults> analyticalResult, std::shared_ptr< ColorConsoleOutput> colorOutput)
{
	return std::shared_ptr<L2NormTest>(new L2NormTest(analyticalResult, colorOutput));
}

void L2NormTest::update()
{
	TestImp::update();
}

void L2NormTest::addSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo)
{
	TestImp::addSimulation(sim, simInfo);
}

void L2NormTest::evaluate()
{
	analyticalResult->calc(simResults.at(0));

	std::vector< double> results = calculator->calc(analyticalResult->getVx(), simResults.at(0)->getVx(), simResults.at(0)->getLevels());

	makeConsoleOutput();
}

std::string L2NormTest::getLogFileOutput()
{
	return std::string();
}

std::vector<bool> L2NormTest::getPassedTests()
{
	return std::vector<bool>();
}

void L2NormTest::makeConsoleOutput()
{

}

L2NormTest::L2NormTest(std::shared_ptr< AnalyticalResults> analyticalResult, std::shared_ptr< ColorConsoleOutput> colorOutput) : TestImp(colorOutput), analyticalResult(analyticalResult)
{
	calculator = L2NormCalculator::getNewInstance();
}