#include "L2NormTestBetweenKernels.h"

#include "Utilities\Calculator\L2NormCalculator\L2NormCalculator.h"
#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"
#include "Utilities\TestSimulation\TestSimulation.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"
#include "Tests\L2NormTestBetweenKernels\PostProcessingStrategy\L2NormBetweenKernelPostProcessingStrategy.h"

#include <iomanip>

std::shared_ptr<L2NormTestBetweenKernels> L2NormTestBetweenKernels::getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep)
{
	return std::shared_ptr<L2NormTestBetweenKernels>(new L2NormTestBetweenKernels(colorOutput, dataToCalculate, timeStep));
}

void L2NormTestBetweenKernels::update()
{
	TestImp::update();
}

void L2NormTestBetweenKernels::evaluate()
{
	basicPostProcessingStrategy->evaluate();
	divergentPostProcessingStrategy->evaluate();

	int tS = calcTimeStepInResults(timeStep);

	if (dataToCalculate == "Vx") {
		basicL2Result = basicPostProcessingStrategy->getL2NormVx(timeStep);
		divergentL2Result = divergentPostProcessingStrategy->getL2NormVx(timeStep);
		resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getVx().at(tS), divergentSimResults->getVx().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getTimeStepLength());
	}	
	if (dataToCalculate == "Vy") {
		basicL2Result = basicPostProcessingStrategy->getL2NormVy(timeStep);
		divergentL2Result = divergentPostProcessingStrategy->getL2NormVy(timeStep);
		resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getVy().at(tS), divergentSimResults->getVy().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getTimeStepLength());
	}
	if (dataToCalculate == "Vz") {
		basicL2Result = basicPostProcessingStrategy->getL2NormVz(timeStep);
		divergentL2Result = divergentPostProcessingStrategy->getL2NormVz(timeStep);
		resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getVz().at(tS), divergentSimResults->getVz().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getTimeStepLength());
	}	
	if (dataToCalculate == "Press") {
		basicL2Result = basicPostProcessingStrategy->getL2NormPress(timeStep);
		divergentL2Result = divergentPostProcessingStrategy->getL2NormPress(timeStep);
		resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getPress().at(tS), divergentSimResults->getPress().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getTimeStepLength());
	}	
	if (dataToCalculate == "Rho") {
		basicL2Result = basicPostProcessingStrategy->getL2NormRho(timeStep);
		divergentL2Result = divergentPostProcessingStrategy->getL2NormRho(timeStep);
		resultL2ToBasicKernel = l2Normcalculator->calc(basicSimResults->getRho().at(tS), divergentSimResults->getRho().at(tS), basicSimResults->getLevels().at(tS), basicSimResults->getNumberOfXNodes(), basicSimResults->getNumberOfZNodes(), basicSimResults->getTimeStepLength());
	}

	if (basicL2Result < divergentL2Result)
			testPassed = true;

	makeConsoleOutput();
}

std::string L2NormTestBetweenKernels::getLogFileOutput()
{
	std::ostringstream oss;
	oss << "L2Norm_" << dataToCalculate << "_TimeStep_" << timeStep << "=" << divergentL2Result << std::endl;
	oss << "L2Norm_BasicKernel_" << dataToCalculate << "_TimeStep_" << timeStep << "=" << basicL2Result << std::endl;
	oss << "L2Norm_Between_Kernels_" << dataToCalculate << "_TimeStep_" << timeStep << "_L" << basicPostProcessingStrategy->getNumberOfXNodes() << "=" << resultL2ToBasicKernel << std::endl << std::endl;
	
	return oss.str();
}

double L2NormTestBetweenKernels::getBasicL2Result()
{
	return basicL2Result;
}

std::vector<bool> L2NormTestBetweenKernels::getPassedTests()
{
	return std::vector<bool>(1, testPassed);
}

void L2NormTestBetweenKernels::makeConsoleOutput()
{
	colorOutput->makeL2NormBetweenKernelsTestOutput(testPassed, basicSimInfo, divergentSimInfo, dataToCalculate, timeStep, basicL2Result, divergentL2Result, resultL2ToBasicKernel);
}

void L2NormTestBetweenKernels::setBasicSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr< L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy)
{
	TestImp::addSimulation(sim, simInfo, postProcessingStrategy);
	this->basicSim = sim;
	this->basicSimInfo = simInfo;
	this->basicPostProcessingStrategy = postProcessingStrategy;
	this->basicSimResults = basicPostProcessingStrategy->getSimulationResult();
}

void L2NormTestBetweenKernels::setDivergentKernelSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy)
{
	TestImp::addSimulation(sim, simInfo, postProcessingStrategy);
	this->divergentSim = sim;
	this->divergentSimInfo = simInfo;
	this->divergentPostProcessingStrategy = postProcessingStrategy;
	this->divergentSimResults = divergentPostProcessingStrategy->getSimulationResult();
}

L2NormTestBetweenKernels::L2NormTestBetweenKernels(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep)
	: TestImp(colorOutput), timeStep(timeStep), dataToCalculate(dataToCalculate)
{
	l2Normcalculator = L2NormCalculator::getInstance();
}

int L2NormTestBetweenKernels::calcTimeStepInResults(unsigned int timeStep)
{
	for (int i = 0; i < basicSimResults->getTimeSteps().size(); i++) {
		if (timeStep == basicSimResults->getTimeSteps().at(i))
			return basicSimResults->getTimeSteps().at(i);
	}
}
