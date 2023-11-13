#include "PostProcessingStrategyL2NormTest.h"

#include "Tests/L2NormTest/L2NormTestParameterStruct.h"

#include "Utilities/Calculator/L2NormCalculator/L2NormCalculatorFactory/L2NormCalculatorFactory.h"
#include "Utilities/Calculator/L2NormCalculator/L2NormCalculator.h"

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"

std::shared_ptr<L2NormPostProcessingStrategy> L2NormPostProcessingStrategy::getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests)
{
    return std::shared_ptr<L2NormPostProcessingStrategy>(new L2NormPostProcessingStrategy(simResult, analyticalResult, testPara, factory, dataToCalcTests));
}

L2NormPostProcessingStrategy::L2NormPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests)
    : PostProcessingStrategyImp(simResult), analyticalResult(analyticalResult), dataToCalculate(dataToCalcTests)
{
    isEvaluated = false;
    basicTimeStep = testPara->basicTimeStep;
    divergentTimeStep = testPara->divergentTimeStep;
    normalizeData = testPara->normalizeData;

    l2NormBasic.resize(dataToCalculate.size());
    l2NormDivergent.resize(dataToCalculate.size());
    for (int i = 0; i < l2NormBasic.size(); i++) {
        l2NormBasic.at(i).resize(normalizeData.size());
        l2NormDivergent.at(i).resize(normalizeData.size());
    }
    
    for (int i = 0; i < normalizeData.size(); i++)
        l2Normcalculator.push_back(factory->makeL2NormCalculator(normalizeData.at(i)));
}

void L2NormPostProcessingStrategy::evaluate()
{
    if (!isEvaluated) {
        analyticalResult->calc(simResult);
        int bS = calcTimeStepInResults(basicTimeStep);
        int dS = calcTimeStepInResults(divergentTimeStep);

        for (int i = 0; i < dataToCalculate.size(); i++) {
            for (int j = 0; j < normalizeData.size(); j++) {
                if (dataToCalculate.at(i) == "Vx") {
                    l2NormBasic.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getVx().at(bS), simResult->getVx().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                    l2NormDivergent.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getVx().at(dS), simResult->getVx().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                }
                if (dataToCalculate.at(i) == "Vy") {
                    l2NormBasic.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getVy().at(bS), simResult->getVy().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                    l2NormDivergent.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getVy().at(dS), simResult->getVy().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                }
                if (dataToCalculate.at(i) == "Vz") {
                    l2NormBasic.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getVz().at(bS), simResult->getVz().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                    l2NormDivergent.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getVz().at(dS), simResult->getVz().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                }
                if (dataToCalculate.at(i) == "Press") {
                    l2NormBasic.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getPress().at(bS), simResult->getPress().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                    l2NormDivergent.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getPress().at(dS), simResult->getPress().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                }
                if (dataToCalculate.at(i) == "Rho") {
                    l2NormBasic.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getRho().at(bS), simResult->getRho().at(bS), simResult->getLevels().at(bS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                    l2NormDivergent.at(i).at(j) = l2Normcalculator.at(j)->calc(analyticalResult->getRho().at(dS), simResult->getRho().at(dS), simResult->getLevels().at(dS), analyticalResult->getNumberOfXNodes(), analyticalResult->getNumberOfZNodes(), analyticalResult->getL0());
                }
            }
        }
        isEvaluated = true;
    }
}

std::vector<double> L2NormPostProcessingStrategy::getL2Norm(std::string aDataToCalc, std::string aNormalizeData)
{
    for (int i = 0; i < dataToCalculate.size(); i++) {
        for (int j = 0; j < normalizeData.size(); j++) {
            if (aDataToCalc == dataToCalculate.at(i) && aNormalizeData == normalizeData.at(j)) {
                std::vector<double> v;
                v.push_back(l2NormBasic.at(i).at(j));
                v.push_back(l2NormDivergent.at(i).at(j));
                return v;
            }
        }
    }

    return std::vector<double>();
}

std::string L2NormPostProcessingStrategy::getErrorMessage(std::string aNormalizeData)
{
    for (int i = 0; i < normalizeData.size(); i++) {
        if (aNormalizeData == normalizeData.at(i))
            return l2Normcalculator.at(i)->getErrorMessage();
    }
    
    return std::string();
}
