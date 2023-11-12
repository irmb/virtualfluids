#ifndef PHIANDNY_TEST_POST_PROCESSING_STRATEGY_H
#define PHIANDNY_TEST_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class FFTCalculator;
struct PhiTestParameterStruct;

class PhiTestPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
    static std::shared_ptr<PhiTestPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests);
    void evaluate();

    double getPhiDiff(std::string dataToCalc);

private:
    PhiTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests);
    
    std::vector<std::vector<double> > reduceDataToTimeSteps(std::vector<std::vector<double> > data);

    std::shared_ptr<AnalyticalResults> analyticalResult;
    std::vector<std::string> dataToCalculate;
    std::shared_ptr<FFTCalculator> fftCalculator;
    unsigned int startTimeStepCalculation;
    unsigned int endTimeStepCalculation;
    std::vector<double> phiDiff;
    bool isEvaluated;
};
#endif 