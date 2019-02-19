#ifndef L2NORM_POST_PROCESSING_STRATEGY_H
#define L2NORM_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class L2NormCalculator;
class L2NormCalculatorFactory;
struct L2NormTestParameterStruct;

class L2NormPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
	static std::shared_ptr<L2NormPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests);
	void evaluate();

	std::vector<double> getL2Norm(std::string dataToCalc, std::string normalizeData);

	std::string getErrorMessage(std::string aNormalizeData);

private:
	L2NormPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara, std::shared_ptr<L2NormCalculatorFactory> factory, std::vector<std::string> dataToCalcTests);
	bool isEvaluated;

	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::vector<std::shared_ptr<L2NormCalculator> > l2Normcalculator;
	
	std::vector<std::string> dataToCalculate;
	std::vector<std::string> normalizeData;
	unsigned int basicTimeStep;
	unsigned int divergentTimeStep;
	std::vector<std::vector<double>> l2NormBasic;
	std::vector<std::vector<double>> l2NormDivergent;
};
#endif