#include "L2NormTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\Calculator\L2NormCalculator\L2NormCalculator.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"

#include <iomanip>

std::shared_ptr<L2NormTest> L2NormTest::getNewInstance(std::shared_ptr< AnalyticalResults> analyticalResult, std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double maxL2NormDiff, unsigned int basicTimeStep, unsigned int divergentTimeStep)
{
	return std::shared_ptr<L2NormTest>(new L2NormTest(analyticalResult, colorOutput, dataToCalculate, maxL2NormDiff, basicTimeStep, divergentTimeStep));
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

	int basicTimeStepInResults = calcTimeStepInResults(basicTimeStep);
	int divergentTimeStepInResults = calcTimeStepInResults(divergentTimeStep);

	resultBasicTimestep = calcL2NormForTimeStep(basicTimeStepInResults);
	resultDivergentTimeStep = calcL2NormForTimeStep(divergentTimeStepInResults);

	diffL2Norm = resultDivergentTimeStep - resultBasicTimestep;

	testPassed = maxL2NormDiff > diffL2Norm;

	makeConsoleOutput();
}

std::string L2NormTest::getLogFileOutput()
{
	std::ostringstream oss;
	oss << "L2Norm_BasicTimeStep_L" << simResults.at(0)->getNumberOfXNodes() << ":" << resultBasicTimestep << std::endl;
	oss << "L2Norm_DivergentTimeStep_L" << simResults.at(0)->getNumberOfXNodes() << ":" << resultDivergentTimeStep << std::endl;
	oss << "L2Norm_Diff_L" << simResults.at(0)->getNumberOfXNodes() << ":" << diffL2Norm << std::endl << std::endl;
	return oss.str();
}

std::vector<bool> L2NormTest::getPassedTests()
{
	return std::vector<bool>(1, testPassed);
}

void L2NormTest::makeConsoleOutput()
{
	colorOutput->makeL2NormTestOutput(testPassed, simInfos.at(0), basicTimeStep, divergentTimeStep, resultBasicTimestep, resultDivergentTimeStep, diffL2Norm);
}

int L2NormTest::calcTimeStepInResults(unsigned int timeStep)
{
	for (int i = 0; i < simResults.at(0)->getTimeSteps().size(); i++) {
		if (timeStep == simResults.at(0)->getTimeSteps().at(i))
			return simResults.at(0)->getTimeSteps().at(i);
	}
}

double L2NormTest::calcL2NormForTimeStep(unsigned int timeStep)
{
	if (dataToCalculate == "Vx")
		return calculator->calc(analyticalResult->getVx().at(timeStep), simResults.at(0)->getVx().at(timeStep), simResults.at(0)->getLevels().at(timeStep));
	if (dataToCalculate == "Vy")
		return calculator->calc(analyticalResult->getVy().at(timeStep), simResults.at(0)->getVy().at(timeStep), simResults.at(0)->getLevels().at(timeStep));
	if (dataToCalculate == "Vz")
		return calculator->calc(analyticalResult->getVz().at(timeStep), simResults.at(0)->getVz().at(timeStep), simResults.at(0)->getLevels().at(timeStep));
	if (dataToCalculate == "Press")
		return calculator->calc(analyticalResult->getPress().at(timeStep), simResults.at(0)->getPress().at(timeStep), simResults.at(0)->getLevels().at(timeStep));
	if (dataToCalculate == "Rho")
		return calculator->calc(analyticalResult->getRho().at(timeStep), simResults.at(0)->getRho().at(timeStep), simResults.at(0)->getLevels().at(timeStep));

}

L2NormTest::L2NormTest(std::shared_ptr< AnalyticalResults> analyticalResult, std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double maxL2NormDiff, unsigned int basicTimeStep, unsigned int divergentTimeStep) : TestImp(colorOutput), analyticalResult(analyticalResult), basicTimeStep(basicTimeStep), divergentTimeStep(divergentTimeStep), dataToCalculate(dataToCalculate), maxL2NormDiff(maxL2NormDiff)
{
	calculator = L2NormCalculator::getNewInstance();
}