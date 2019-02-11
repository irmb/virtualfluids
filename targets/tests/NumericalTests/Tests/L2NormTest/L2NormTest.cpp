#include "L2NormTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"

#include "Tests/L2NormTest/PostProcessingStrategy/PostProcessingStrategyL2NormTest.h"
#include "Tests/L2NormTest/L2NormTestParameterStruct.h"

#include <iomanip>

std::shared_ptr<L2NormTest> L2NormTest::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData)
{
	return std::shared_ptr<L2NormTest>(new L2NormTest(colorOutput, testParameter, dataToCalculate, maxL2NormDiff, normalizeData));
}

void L2NormTest::update()
{
	TestImp::update();
}

void L2NormTest::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormPostProcessingStrategy> postProStrategy)
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

	if (resultBasicTimestep < 0 || resultDivergentTimeStep < 0) {
		testError = true;
		testPassed = false;
	}
	else
	{
		testPassed = maxL2NormDiff > diffL2Norm;
	}
	

	makeConsoleOutput();
}

std::string L2NormTest::getLogFileOutput()
{
	std::ostringstream oss;
	oss << "NormalizeData_L" << l2NormPostProStrategies.at(0)->getNumberOfXNodes() << "_" << dataToCalculate << "=" << normalizeData << std::endl;
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
	if (!testError)
		colorOutput->makeL2NormTestOutput(testPassed, simInfos.at(0), normalizeData, basicTimeStep, divergentTimeStep, dataToCalculate, resultBasicTimestep, resultDivergentTimeStep, diffL2Norm);
	else
		colorOutput->makeL2NormTestErrorOutput(l2NormPostProStrategies.at(0)->getErrorMessage(), simInfos.at(0), normalizeData, basicTimeStep, divergentTimeStep, dataToCalculate);

}

L2NormTest::L2NormTest(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData)
	: TestImp(colorOutput), dataToCalculate(dataToCalculate), normalizeData(normalizeData)
{
	basicTimeStep = testParameter->basicTimeStep;
	divergentTimeStep = testParameter->divergentTimeStep;
	this->maxL2NormDiff = maxL2NormDiff;
	testError = false;
}