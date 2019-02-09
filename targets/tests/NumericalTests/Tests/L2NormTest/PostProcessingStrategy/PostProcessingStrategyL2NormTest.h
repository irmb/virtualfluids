#ifndef L2NORM_POST_PROCESSING_STRATEGY_H
#define L2NORM_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class L2NormCalculator;
struct L2NormTestParameterStruct;

class L2NormPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
	static std::shared_ptr<L2NormPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara);
	void evaluate();

	std::vector<double> getL2NormVx();
	std::vector<double> getL2NormVy();
	std::vector<double> getL2NormVz();
	std::vector<double> getL2NormPress();
	std::vector<double> getL2NormRho();

private:
	L2NormPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestParameterStruct> testPara);
	bool isEvaluated;

	std::shared_ptr<L2NormCalculator> l2Normcalculator;
	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::vector<std::string> dataToCalculateL2;
	unsigned int basicTimeStepL2Norm;
	unsigned int divergentTimeStepL2Norm;
	double l2VxBasic, l2VyBasic, l2VzBasic, l2RhoBasic, l2PressBasic;
	double l2VxDivergent, l2VyDivergent, l2VzDivergent, l2RhoDivergent, l2PressDivergent;
};
#endif