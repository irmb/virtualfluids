#ifndef NY_TEST_POST_PROCESSING_STRATEGY_H
#define NY_TEST_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class FFTCalculator;
struct NyTestParameterStruct;

class NyTestPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
    static std::shared_ptr<NyTestPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<NyTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests);
    void evaluate();

    double getNy(std::string dataToCalculate);

private:
    NyTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<NyTestParameterStruct> testPara, std::vector<std::string> dataToCalcTests);
    
    std::vector<std::vector<double> > reduceDataToTimeSteps(std::vector<std::vector<double> > data, unsigned int startTimeStep, unsigned int endTimeStep);

    std::shared_ptr<AnalyticalResults> analyticalResult;
    std::vector<std::string> dataToCalculate;
    std::shared_ptr<FFTCalculator> fftCalculator;
    unsigned int startTimeStepCalculation;
    unsigned int endTimeStepCalculation;
    std::vector<double> ny;
    bool isEvaluated;
};
#endif 