#include "L2NormBetweenKernelPostProcessingStrategy.h"

#include "Utilities/Calculator/L2NormCalculator/L2NormCalculator.h"
#include "Utilities/Calculator/L2NormCalculator/L2NormCalculatorFactory/L2NormCalculatorFactory.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"

#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsParameterStruct.h"

std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> L2NormBetweenKernelPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests)
{
	return std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy>(new L2NormBetweenKernelPostProcessingStrategy(simResult, analyticalResult, testPara, factory, dataToCalcTests));
}

void L2NormBetweenKernelPostProcessingStrategy::evaluate()
{
	if (!isEvaluated) {
		analyticalResult->calc(simResult);

		l2Norm.resize(dataToCalculate.size());
		for (int i = 0; i < dataToCalculate.size(); i++) {
			l2Norm.at(i).resize(normalizeData.size());
			for (int j = 0; j < normalizeData.size(); j++) {
				l2Norm.at(i).at(j).resize(timeSteps.size());
			}
		}

		for (int i = 0; i < dataToCalculate.size(); i++) {
			for (int j = 0; j < normalizeData.size(); j++) {
				for (int k = 0; k < timeSteps.size(); k++) {
					int time = calcTimeStepInResults(timeSteps.at(k));
					if (dataToCalculate.at(i) == "Vx")
						l2Norm.at(i).at(j).at(k) = l2Normcalculator.at(j)->calc(analyticalResult->getVx().at(time), simResult->getVx().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes());
					if (dataToCalculate.at(i) == "Vy")
						l2Norm.at(i).at(j).at(k) = l2Normcalculator.at(j)->calc(analyticalResult->getVy().at(time), simResult->getVy().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes());
					if (dataToCalculate.at(i) == "Vz")
						l2Norm.at(i).at(j).at(k) = l2Normcalculator.at(j)->calc(analyticalResult->getVz().at(time), simResult->getVz().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes());
					if (dataToCalculate.at(i) == "Press")
						l2Norm.at(i).at(j).at(k) = l2Normcalculator.at(j)->calc(analyticalResult->getPress().at(time), simResult->getPress().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes());
					if (dataToCalculate.at(i) == "Rho")
						l2Norm.at(i).at(j).at(k) = l2Normcalculator.at(j)->calc(analyticalResult->getRho().at(time), simResult->getRho().at(time), simResult->getLevels().at(time), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes());
				}
			}
		}
		isEvaluated = true;
	}	
}

double L2NormBetweenKernelPostProcessingStrategy::getL2Norm(std::string aDataToCalc, std::string aNormalizeData, int aTimeStep)
{
	for (int i = 0; i < dataToCalculate.size(); i++) {
		for (int j = 0; j < normalizeData.size(); j++) {
			for (int k = 0; k < timeSteps.size(); k++) {
				if (aDataToCalc == dataToCalculate.at(i) && aNormalizeData == normalizeData.at(j) && aTimeStep == timeSteps.at(k))
					return l2Norm.at(i).at(j).at(k);
			}
		}
	}

	return 0.0;
}

std::string L2NormBetweenKernelPostProcessingStrategy::getErrorMessage(std::string aNormalizeData)
{
	for (int i = 0; i < normalizeData.size(); i++) {
		if (aNormalizeData == normalizeData.at(i))
			return l2Normcalculator.at(i)->getErrorMessage();
	}
	return std::string();
}

std::shared_ptr<SimulationResults> L2NormBetweenKernelPostProcessingStrategy::getSimulationResult()
{
	return simResult;
}

L2NormBetweenKernelPostProcessingStrategy::L2NormBetweenKernelPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests)
	: PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult), dataToCalculate(dataToCalcTests)
{
	isEvaluated = false;
	normalizeData = testPara->normalizeData;
	timeSteps = testPara->timeSteps;

	l2Norm.resize(dataToCalculate.size());
	for (int i = 0; i < dataToCalculate.size(); i++) {
		l2Norm.at(i).resize(normalizeData.size());
		for (int j = 0; j < normalizeData.size(); j++) {
			l2Norm.at(i).at(j).resize(timeSteps.size());
		}
	}


	for (int i = 0; i < normalizeData.size(); i++)
		l2Normcalculator.push_back(factory->makeL2NormCalculator(normalizeData.at(i)));
	
}

int L2NormBetweenKernelPostProcessingStrategy::calcPosInTimeStep(int time)
{
	for (int i = 0; i < timeSteps.size(); i++)
		if (timeSteps.at(i) == time)
			return i;
}
