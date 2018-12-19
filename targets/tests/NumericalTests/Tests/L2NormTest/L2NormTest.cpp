#include "L2NormTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\PostProcessingResults\PostProcessingResults.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"

#include "Tests\L2NormTest\PostProcessingStrategy\PostProcessingStrategyL2NormTest.h"

#include <iomanip>

std::shared_ptr<L2NormTest> L2NormTest::getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double maxL2NormDiff, unsigned int basicTimeStep, unsigned int divergentTimeStep)
{
	return std::shared_ptr<L2NormTest>(new L2NormTest(colorOutput, dataToCalculate, maxL2NormDiff, basicTimeStep, divergentTimeStep));
}

void L2NormTest::update()
{
	TestImp::update();
}

void L2NormTest::addSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr< L2NormPostProcessingStrategy> postProStrategy)
{
	TestImp::addSimulation(sim, simInfo, postProStrategy);
	l2NormPostProStrategies.push_back(postProStrategy);
}

void L2NormTest::evaluate()
{
	std::vector<double> results;

	if (dataToCalculate == "Vx")
		results = l2NormPostProStrategies.at(0)->getL2NormVx();
	if (dataToCalculate == "Vy")
		results = l2NormPostProStrategies.at(0)->getL2NormVy();
	if (dataToCalculate == "Vz")
		results = l2NormPostProStrategies.at(0)->getL2NormVz();
	if (dataToCalculate == "Press")
		results = l2NormPostProStrategies.at(0)->getL2NormPress();
	if (dataToCalculate == "Rho")
		results = l2NormPostProStrategies.at(0)->getL2NormRho();
		
	resultBasicTimestep = results.at(0);
	resultDivergentTimeStep = results.at(1);
	diffL2Norm = resultDivergentTimeStep - resultBasicTimestep;

	testPassed = maxL2NormDiff > diffL2Norm;

	makeConsoleOutput();
}

std::string L2NormTest::getLogFileOutput()
{
	std::ostringstream oss;
	oss << "L2Norm_BasicTimeStep_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "=" << resultBasicTimestep << std::endl;
	oss << "L2Norm_DivergentTimeStep_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "=" << resultDivergentTimeStep << std::endl;
	oss << "L2Norm_Diff_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "=" << diffL2Norm << std::endl << std::endl;
	return oss.str();
}

std::vector<bool> L2NormTest::getPassedTests()
{
	return std::vector<bool>(1, testPassed);
}

void L2NormTest::makeConsoleOutput()
{
	colorOutput->makeL2NormTestOutput(testPassed, simInfos.at(0), basicTimeStep, divergentTimeStep, dataToCalculate, resultBasicTimestep, resultDivergentTimeStep, diffL2Norm);
}

L2NormTest::L2NormTest(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double maxL2NormDiff, unsigned int basicTimeStep, unsigned int divergentTimeStep)
	: TestImp(colorOutput), basicTimeStep(basicTimeStep), divergentTimeStep(divergentTimeStep), dataToCalculate(dataToCalculate), maxL2NormDiff(maxL2NormDiff)
{
	
}