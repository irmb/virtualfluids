#include "PostProcessingStrategyPhiAndNuTest.h"

#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Calculator\FFTCalculator\FFTCalculator.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"

#include "Tests\PhiAndNuTest\PhiAndNuTestParameterStruct.h"

std::shared_ptr<PhiAndNuTestPostProcessingStrategy> PhiAndNuTestPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiAndNuTestParameterStruct> testPara)
{
	return std::shared_ptr<PhiAndNuTestPostProcessingStrategy>(new PhiAndNuTestPostProcessingStrategy(simResult, analyticalResult, testPara));
}

PhiAndNuTestPostProcessingStrategy::PhiAndNuTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiAndNuTestParameterStruct> testPara)
	: PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult)
{
	dataToCalculatePhiAndNu = testPara->basicTestParameter->dataToCalc;
	startTimeStepCalculationPhiNu = testPara->startTimeStepCalculation;
	endTimeStepCalculationPhiNu = testPara->endTimeStepCalculation;

	isEvaluated = false;
	fftCalculator = FFTCalculator::getNewInstance(simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
}

std::vector<std::vector<double>> PhiAndNuTestPostProcessingStrategy::reduceDataToTimeSteps(std::vector<std::vector<double>> data, unsigned int startTimeStep, unsigned int endTimeStep)
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

void PhiAndNuTestPostProcessingStrategy::evaluate()
{
	if (!isEvaluated) {
		analyticalResult->calc(simResult);
		for (int i = 0; i < dataToCalculatePhiAndNu.size(); i++) {
			if (dataToCalculatePhiAndNu.at(i) == "Vx") {
				fftCalculator->calc(reduceDataToTimeSteps(simResult->getVx(), startTimeStepCalculationPhiNu, endTimeStepCalculationPhiNu), false);
				nuVx = fftCalculator->getNu();
				phiDiffVx = fftCalculator->getPhiDiff();
			}
			if (dataToCalculatePhiAndNu.at(i) == "Vz") {
				fftCalculator->calc(reduceDataToTimeSteps(simResult->getVz(), startTimeStepCalculationPhiNu, endTimeStepCalculationPhiNu), true);
				nuVz = fftCalculator->getNu();
				phiDiffVz = fftCalculator->getPhiDiff();
			}
		}
	}
}

double PhiAndNuTestPostProcessingStrategy::getNuVx()
{
	return nuVx;
}

double PhiAndNuTestPostProcessingStrategy::getNuVy()
{
	return nuVy;
}

double PhiAndNuTestPostProcessingStrategy::getNuVz()
{
	return nuVz;
}

double PhiAndNuTestPostProcessingStrategy::getPhiDiffVx()
{
	return phiDiffVx;
}

double PhiAndNuTestPostProcessingStrategy::getPhiDiffVy()
{
	return phiDiffVy;
}

double PhiAndNuTestPostProcessingStrategy::getPhiDiffVz()
{
	return phiDiffVz;
}