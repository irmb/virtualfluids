#include "PostProcessingStrategyL2NormTest.h"

#include "Tests\L2NormTest\L2NormTestParameterStruct.h"

#include "Utilities\Calculator\L2NormCalculator\L2NormCalculator.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Results\SimulationResults\SimulationResults.h"

std::shared_ptr<L2NormPostProcessingStrategy> L2NormPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara)
{
	return std::shared_ptr<L2NormPostProcessingStrategy>(new L2NormPostProcessingStrategy(simResult, analyticalResult, testPara));
}

L2NormPostProcessingStrategy::L2NormPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara)
	: PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult)
{
	dataToCalculateL2 = testPara->basicTestParameter->dataToCalc;
	basicTimeStepL2Norm = testPara->basicTimeStep;
	divergentTimeStepL2Norm = testPara->divergentTimeStep;
	
	isEvaluated = false;
	l2Normcalculator = L2NormCalculator::getInstance();
}

void L2NormPostProcessingStrategy::evaluate()
{
	if (!isEvaluated) {
		analyticalResult->calc(simResult);
		int bS = calcTimeStepInResults(basicTimeStepL2Norm);
		int dS = calcTimeStepInResults(divergentTimeStepL2Norm);

		for (int i = 0; i < dataToCalculateL2.size(); i++) {
			if (dataToCalculateL2.at(i) == "Vx") {
				l2VxBasic = l2Normcalculator->calc(analyticalResult->getVx().at(bS), simResult->getVx().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
				l2VxDivergent = l2Normcalculator->calc(analyticalResult->getVx().at(dS), simResult->getVx().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
			}
			if (dataToCalculateL2.at(i) == "Vy") {
				l2VyBasic = l2Normcalculator->calc(analyticalResult->getVy().at(bS), simResult->getVy().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
				l2VyDivergent = l2Normcalculator->calc(analyticalResult->getVy().at(dS), simResult->getVy().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
			}
			if (dataToCalculateL2.at(i) == "Vz") {
				l2VzBasic = l2Normcalculator->calc(analyticalResult->getVz().at(bS), simResult->getVz().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
				l2VzDivergent = l2Normcalculator->calc(analyticalResult->getVz().at(dS), simResult->getVz().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
			}
			if (dataToCalculateL2.at(i) == "Press") {
				l2PressBasic = l2Normcalculator->calc(analyticalResult->getPress().at(bS), simResult->getPress().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
				l2PressDivergent = l2Normcalculator->calc(analyticalResult->getPress().at(dS), simResult->getPress().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
			}
			if (dataToCalculateL2.at(i) == "Rho") {
				l2RhoBasic = l2Normcalculator->calc(analyticalResult->getRho().at(bS), simResult->getRho().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
				l2RhoDivergent = l2Normcalculator->calc(analyticalResult->getRho().at(dS), simResult->getRho().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getTimeStepLength());
			}
		}
		isEvaluated = true;
	}
}

std::vector<double> L2NormPostProcessingStrategy::getL2NormVx()
{
	std::vector<double> v;
	v.push_back(l2VxBasic);
	v.push_back(l2VxDivergent);
	return v;
}

std::vector<double> L2NormPostProcessingStrategy::getL2NormVy()
{
	std::vector<double> v;
	v.push_back(l2VyBasic);
	v.push_back(l2VyDivergent);
	return v;
}

std::vector<double> L2NormPostProcessingStrategy::getL2NormVz()
{
	std::vector<double> v;
	v.push_back(l2VzBasic);
	v.push_back(l2VzDivergent);
	return v;
}

std::vector<double> L2NormPostProcessingStrategy::getL2NormPress()
{
	std::vector<double> v;
	v.push_back(l2PressBasic);
	v.push_back(l2PressDivergent);
	return v;
}

std::vector<double> L2NormPostProcessingStrategy::getL2NormRho()
{
	std::vector<double> v;
	v.push_back(l2RhoBasic);
	v.push_back(l2RhoDivergent);
	return v;
}