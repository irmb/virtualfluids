#include "PostProcessingResultsImp.h"

#include "Utilities\Calculator\FFTCalculator\FFTCalculator.h"
#include "Utilities\Calculator\L2NormCalculator\L2NormCalculator.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"

std::shared_ptr<PostProcessingResultsImp> PostProcessingResultsImp::getNewInstance(std::shared_ptr< SimulationResults> simResult, std::shared_ptr< AnalyticalResults> analyticalResult, std::vector<std::string> dataToCalculatePhiAndNu, std::vector<std::string> dataToCalculateL2, unsigned int startTimeStepCalculationPhiNu, unsigned int endTimeStepCalculationPhiNu, unsigned int basicTimeStepL2Norm, unsigned int divergentTimeStepL2Norm)
{
	return std::shared_ptr<PostProcessingResultsImp>(new PostProcessingResultsImp(simResult, analyticalResult, dataToCalculatePhiAndNu, dataToCalculateL2, startTimeStepCalculationPhiNu, endTimeStepCalculationPhiNu, basicTimeStepL2Norm, divergentTimeStepL2Norm));
}

PostProcessingResultsImp::PostProcessingResultsImp(std::shared_ptr< SimulationResults> simResult, std::shared_ptr< AnalyticalResults> analyticalResult, std::vector<std::string> dataToCalculatePhiAndNu, std::vector<std::string> dataToCalculateL2, unsigned int startTimeStepCalculationPhiNu, unsigned int endTimeStepCalculationPhiNu, unsigned int basicTimeStepL2Norm, unsigned int divergentTimeStepL2Norm)
	: simResult(simResult), analyticalResult(analyticalResult), dataToCalculatePhiAndNu(dataToCalculatePhiAndNu), dataToCalculateL2(dataToCalculateL2), startTimeStepCalculationPhiNu(startTimeStepCalculationPhiNu), endTimeStepCalculationPhiNu(endTimeStepCalculationPhiNu), basicTimeStepL2Norm(basicTimeStepL2Norm), divergentTimeStepL2Norm(divergentTimeStepL2Norm)
{
	isEvaluated = false;
	fftCalculator = FFTCalculator::getNewInstance(simResult->getNumberOfXNodes(), simResult->getNumberOfZNodes(), simResult->getTimeStepLength());
	l2Normcalculator = L2NormCalculator::getNewInstance();
}

void PostProcessingResultsImp::evaluate()
{
	if (!isEvaluated) {
		analyticalResult->calc(simResult);
		for (int i = 0; i < dataToCalculatePhiAndNu.size(); i++) {
			if (dataToCalculatePhiAndNu.at(i) == "Vx") {
				fftCalculator->calc(reduceDataToTimeSteps(simResult->getVx(), startTimeStepCalculationPhiNu, endTimeStepCalculationPhiNu));
				nuVx = fftCalculator->getNu();
				phiDiffVx = fftCalculator->getPhiDiff();
			}
			if (dataToCalculatePhiAndNu.at(i) == "Vz") {
				fftCalculator->calc(reduceDataToTimeSteps(simResult->getVz(), startTimeStepCalculationPhiNu, endTimeStepCalculationPhiNu));
				nuVz = fftCalculator->getNu();
				phiDiffVz = fftCalculator->getPhiDiff();
			}
		}
		analyticalResult->calc(simResult);
		int bS = calcTimeStepInResults(basicTimeStepL2Norm);
		int dS = calcTimeStepInResults(divergentTimeStepL2Norm);

		for (int i = 0; i < dataToCalculateL2.size(); i++) {
			if (dataToCalculateL2.at(i) == "Vx"){
				l2VxBasic = l2Normcalculator->calc(analyticalResult->getVx().at(bS), simResult->getVx().at(bS), simResult->getLevels().at(bS));
				l2VxDivergent = l2Normcalculator->calc(analyticalResult->getVx().at(dS), simResult->getVx().at(dS), simResult->getLevels().at(dS));
			}
			if (dataToCalculateL2.at(i) == "Vy") {
				l2VyBasic = l2Normcalculator->calc(analyticalResult->getVy().at(bS), simResult->getVy().at(bS), simResult->getLevels().at(bS));
				l2VyDivergent = l2Normcalculator->calc(analyticalResult->getVy().at(dS), simResult->getVy().at(dS), simResult->getLevels().at(dS));
			}
			if (dataToCalculateL2.at(i) == "Vz") {
				l2VzBasic = l2Normcalculator->calc(analyticalResult->getVz().at(bS), simResult->getVz().at(bS), simResult->getLevels().at(bS));
				l2VzDivergent = l2Normcalculator->calc(analyticalResult->getVz().at(dS), simResult->getVz().at(dS), simResult->getLevels().at(dS));
			}
			if (dataToCalculateL2.at(i) == "Press") {
				l2PressBasic = l2Normcalculator->calc(analyticalResult->getPress().at(bS), simResult->getPress().at(bS), simResult->getLevels().at(bS));
				l2PressDivergent = l2Normcalculator->calc(analyticalResult->getPress().at(dS), simResult->getPress().at(dS), simResult->getLevels().at(dS));
			}
			if (dataToCalculateL2.at(i) == "Rho") {
				l2RhoBasic = l2Normcalculator->calc(analyticalResult->getRho().at(bS), simResult->getRho().at(bS), simResult->getLevels().at(bS));
				l2RhoDivergent = l2Normcalculator->calc(analyticalResult->getRho().at(dS), simResult->getRho().at(dS), simResult->getLevels().at(dS));
			}
		}
		isEvaluated = true;
	}
}

double PostProcessingResultsImp::getNuVx()
{
	return nuVx;
}

double PostProcessingResultsImp::getNuVy()
{
	return nuVy;
}

double PostProcessingResultsImp::getNuVz()
{
	return nuVz;
}

double PostProcessingResultsImp::getPhiDiffVx()
{
	return phiDiffVx;
}

double PostProcessingResultsImp::getPhiDiffVy()
{
	return phiDiffVy;
}

double PostProcessingResultsImp::getPhiDiffVz()
{
	return phiDiffVz;
}

std::vector<double> PostProcessingResultsImp::getL2NormVx()
{
	std::vector<double> v;
	v.push_back(l2VxBasic);
	v.push_back(l2VxDivergent);
	return v;
}

std::vector<double> PostProcessingResultsImp::getL2NormVy()
{
	std::vector<double> v;
	v.push_back(l2VyBasic);
	v.push_back(l2VyDivergent);
	return v;
}

std::vector<double> PostProcessingResultsImp::getL2NormVz()
{
	std::vector<double> v;
	v.push_back(l2VzBasic);
	v.push_back(l2VzDivergent);
	return v;
}

std::vector<double> PostProcessingResultsImp::getL2NormPress()
{
	std::vector<double> v;
	v.push_back(l2PressBasic);
	v.push_back(l2PressDivergent);
	return v;
}

std::vector<double> PostProcessingResultsImp::getL2NormRho()
{
	std::vector<double> v;
	v.push_back(l2RhoBasic);
	v.push_back(l2RhoDivergent);
	return v;
}

int PostProcessingResultsImp::getNumberOfXNodes()
{
	return simResult->getNumberOfXNodes();
}

int PostProcessingResultsImp::getNumberOfYNodes()
{
	return simResult->getNumberOfYNodes();
}

int PostProcessingResultsImp::getNumberOfZNodes()
{
	return simResult->getNumberOfZNodes();
}

std::vector<std::vector<double>> PostProcessingResultsImp::reduceDataToTimeSteps(std::vector<std::vector<double>> data, unsigned int startTimeStep, unsigned int endTimeStep)
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

int PostProcessingResultsImp::calcTimeStepInResults(unsigned int timeStep)
{
	for (int i = 0; i < simResult->getTimeSteps().size(); i++) {
		if (timeStep == simResult->getTimeSteps().at(i))
			return simResult->getTimeSteps().at(i);
	}
}
