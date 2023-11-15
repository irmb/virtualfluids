#ifndef L2NORM_BETWEEN_KERNEL_POST_PROCESSING_STRATEGY_H
#define L2NORM_BETWEEN_KERNEL_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class L2NormCalculator;
class L2NormCalculatorFactory;
struct L2NormTestBetweenKernelsParameterStruct;

class L2NormBetweenKernelPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
    static std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests);
    void evaluate();

    double getL2Norm(std::string aDataToCalc, std::string aNormalizeData, int aTimeStep);

    std::string getErrorMessage(std::string aNormalizeData);

    virtual std::shared_ptr<SimulationResults> getSimulationResult();

private:
    L2NormBetweenKernelPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests);

    int calcPosInTimeStep(int time);

    std::shared_ptr<AnalyticalResults> analyticalResult;
    std::vector<std::shared_ptr<L2NormCalculator> > l2Normcalculator;
    std::vector<std::string> dataToCalculate;
    std::vector<std::string> normalizeData;
    std::vector<int> timeSteps;

    std::vector<std::vector<std::vector<double> > > l2Norm;    
    bool isEvaluated;

};
#endif