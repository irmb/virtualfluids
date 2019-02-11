#include "L2NormBetweenKernelPostProcessingStrategy.h"

#include "Utilities/Calculator/L2NormCalculator/L2NormCalculator.h"
#include "Utilities/Calculator/L2NormCalculator/L2NormCalculatorFactory/L2NormCalculatorFactory.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"

#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsParameterStruct.h"

std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> L2NormBetweenKernelPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara)
{
	return std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy>(new L2NormBetweenKernelPostProcessingStrategy(simResult, analyticalResult, testPara));
}

void L2NormBetweenKernelPostProcessingStrategy::evaluate()
{
	if (!isEvaluated) {
		analyticalResult->calc(simResult);

		l2Norm.resize(timeSteps.size());
		for (int j = 0; j < dataToCalculate.size(); j++) {
			for (int i = 0; i < timeSteps.size(); i++) {
				int time = calcTimeStepInResults(timeSteps.at(i));
				if (dataToCalculate.at(j) == "Vx")
					l2Vx.push_back(l2Normcalculator->calc(analyticalResult->getVx().at(time), simResult->getVx().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength()));
				if (dataToCalculate.at(j) == "Vy")
					l2Vy.push_back(l2Normcalculator->calc(analyticalResult->getVy().at(time), simResult->getVy().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength()));
				if (dataToCalculate.at(j) == "Vz")
					l2Vz.push_back(l2Normcalculator->calc(analyticalResult->getVz().at(time), simResult->getVz().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength()));
				if (dataToCalculate.at(j) == "Press")
					l2Press.push_back(l2Normcalculator->calc(analyticalResult->getPress().at(time), simResult->getPress().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength()));
				if (dataToCalculate.at(j) == "Rho")
					l2Rho.push_back(l2Normcalculator->calc(analyticalResult->getRho().at(time), simResult->getRho().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength()));
			}
		}
		isEvaluated = true;
	}	
}

double L2NormBetweenKernelPostProcessingStrategy::getL2NormVx(int timeStep)
{
	return l2Vx.at(calcPosInTimeStep(timeStep));
}

double L2NormBetweenKernelPostProcessingStrategy::getL2NormVy(int timeStep)
{
	return l2Vy.at(calcPosInTimeStep(timeStep));
}

double L2NormBetweenKernelPostProcessingStrategy::getL2NormVz(int timeStep)
{
	return l2Vz.at(calcPosInTimeStep(timeStep));
}

double L2NormBetweenKernelPostProcessingStrategy::getL2NormPress(int timeStep)
{
	return l2Press.at(calcPosInTimeStep(timeStep));
}

double L2NormBetweenKernelPostProcessingStrategy::getL2NormRho(int timeStep)
{
	return l2Rho.at(calcPosInTimeStep(timeStep));
}

std::string L2NormBetweenKernelPostProcessingStrategy::getErrorMessage()
{
	return l2Normcalculator->getErrorMessage();
}

std::shared_ptr<SimulationResults> L2NormBetweenKernelPostProcessingStrategy::getSimulationResult()
{
	return simResult;
}

L2NormBetweenKernelPostProcessingStrategy::L2NormBetweenKernelPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara)
	: PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult)
{
	timeSteps = testPara->timeSteps;
	dataToCalculate = testPara->basicTestParameter->dataToCalc;
	std::shared_ptr<L2NormCalculatorFactory> l2NormCalculatorFactory = L2NormCalculatorFactory::getInstance();
	l2Normcalculator = l2NormCalculatorFactory->makeL2NormCalculator(testPara->normalizeWith);
	isEvaluated = false;
}

int L2NormBetweenKernelPostProcessingStrategy::calcPosInTimeStep(int time)
{
	for (int i = 0; i < timeSteps.size(); i++)
		if (timeSteps.at(i) == time)
			return i;
}
