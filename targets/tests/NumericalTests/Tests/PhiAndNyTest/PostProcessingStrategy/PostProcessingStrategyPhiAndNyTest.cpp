#include "PostProcessingStrategyPhiAndNyTest.h"

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"

#include "Tests/PhiAndNyTest/PhiAndNyTestParameterStruct.h"

std::shared_ptr<PhiAndNyTestPostProcessingStrategy> PhiAndNyTestPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiAndNyTestParameterStruct> testPara)
{
	return std::shared_ptr<PhiAndNyTestPostProcessingStrategy>(new PhiAndNyTestPostProcessingStrategy(simResult, analyticalResult, testPara));
}

PhiAndNyTestPostProcessingStrategy::PhiAndNyTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiAndNyTestParameterStruct> testPara)
	: PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult)
{
	dataToCalculatePhiAndNy = testPara->basicTestParameter->dataToCalc;
	startTimeStepCalculationPhiNy = testPara->startTimeStepCalculation;
	endTimeStepCalculationPhiNy = testPara->endTimeStepCalculation;

	isEvaluated = false;
	fftCalculator = FFTCalculator::getNewInstance(simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
}

std::vector<std::vector<double> > PhiAndNyTestPostProcessingStrategy::reduceDataToTimeSteps(std::vector<std::vector<double> > data, unsigned int startTimeStep, unsigned int endTimeStep)
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

void PhiAndNyTestPostProcessingStrategy::evaluate()
{
	if (!isEvaluated) {
		analyticalResult->calc(simResult);
		for (int i = 0; i < dataToCalculatePhiAndNy.size(); i++) {
			if (dataToCalculatePhiAndNy.at(i) == "Vx") {
				fftCalculator->calc(reduceDataToTimeSteps(simResult->getVx(), startTimeStepCalculationPhiNy, endTimeStepCalculationPhiNy), false);
				nyVx = fftCalculator->getNy();
				phiDiffVx = fftCalculator->getPhiDiff();
			}
			if (dataToCalculatePhiAndNy.at(i) == "Vz") {
				fftCalculator->calc(reduceDataToTimeSteps(simResult->getVz(), startTimeStepCalculationPhiNy, endTimeStepCalculationPhiNy), true);
				nyVz = fftCalculator->getNy();
				phiDiffVz = fftCalculator->getPhiDiff();
			}
		}
	}
}

double PhiAndNyTestPostProcessingStrategy::getNyVx()
{
	return nyVx;
}

double PhiAndNyTestPostProcessingStrategy::getNyVy()
{
	return nyVy;
}

double PhiAndNyTestPostProcessingStrategy::getNyVz()
{
	return nyVz;
}

double PhiAndNyTestPostProcessingStrategy::getPhiDiffVx()
{
	return phiDiffVx;
}

double PhiAndNyTestPostProcessingStrategy::getPhiDiffVy()
{
	return phiDiffVy;
}

double PhiAndNyTestPostProcessingStrategy::getPhiDiffVz()
{
	return phiDiffVz;
}