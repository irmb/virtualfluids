#include "PhiTestPostProcessingStrategy.h"

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"

#include "Tests/PhiTest/PhiTestParameterStruct.h"

std::shared_ptr<PhiTestPostProcessingStrategy> PhiTestPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
{
	return std::shared_ptr<PhiTestPostProcessingStrategy>(new PhiTestPostProcessingStrategy(simResult, analyticalResult, testPara, dataToCalcTests));
}

PhiTestPostProcessingStrategy::PhiTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
	: PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult), dataToCalculate(dataToCalcTests)
{
	startTimeStepCalculation = testPara->startTimeStepCalculation;
	endTimeStepCalculation = testPara->endTimeStepCalculation;
	phiDiff.resize(dataToCalculate.size());

	isEvaluated = false;
	fftCalculator = FFTCalculator::getInstance();
}

std::vector<std::vector<double> > PhiTestPostProcessingStrategy::reduceDataToTimeSteps(std::vector<std::vector<double> > data, unsigned int startTimeStep, unsigned int endTimeStep)
{
	std::vector<int> timeStepsToDelete;

	for (int i = simResult->getTimeSteps().size() - 1; i >= 0; i--) {
		if (simResult->getTimeSteps().at(i) > endTimeStep)
			timeStepsToDelete.push_back(i);
		if (simResult->getTimeSteps().at(i) < startTimeStep)
			timeStepsToDelete.push_back(i);
	}

	for (int i = 0; i < timeStepsToDelete.size(); i++)
		data.erase(data.begin() + timeStepsToDelete.at(i));

	return data;
}

void PhiTestPostProcessingStrategy::evaluate()
{
	if (!isEvaluated) {
		for (int i = 0; i < dataToCalculate.size(); i++) {
			if (dataToCalculate.at(i) == "Vx") {
				phiDiff.at(i) = fftCalculator->calcPhiDiff(reduceDataToTimeSteps(simResult->getVx(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
			}
			if (dataToCalculate.at(i) == "Vy") {
				phiDiff.at(i) = abs(fftCalculator->calcPhiDiff(reduceDataToTimeSteps(simResult->getVy(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength()));
			}
			if (dataToCalculate.at(i) == "Vz") {
				phiDiff.at(i) = fftCalculator->calcPhiDiff(reduceDataToTimeSteps(simResult->getVz(), startTimeStepCalculation, endTimeStepCalculation), true, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
			}
			if (dataToCalculate.at(i) == "Press") {
				phiDiff.at(i) = fftCalculator->calcPhiDiff(reduceDataToTimeSteps(simResult->getPress(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
			}
			if (dataToCalculate.at(i) == "Rho") {
				phiDiff.at(i) = fftCalculator->calcPhiDiff(reduceDataToTimeSteps(simResult->getRho(), startTimeStepCalculation, endTimeStepCalculation), false, simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
			}
		}
		isEvaluated = true;
	}
}

double PhiTestPostProcessingStrategy::getPhiDiff(std::string dataToCalc)
{
	for (int i = 0; i < dataToCalculate.size(); i++)
		if(dataToCalculate.at(i) == dataToCalc)
			return phiDiff.at(i);
}
