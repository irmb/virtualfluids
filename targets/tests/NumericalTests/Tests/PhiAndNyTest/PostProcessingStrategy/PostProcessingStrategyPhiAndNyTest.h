#ifndef PHIANDNY_TEST_POST_PROCESSING_STRATEGY_H
#define PHIANDNY_TEST_POST_PROCESSING_STRATEGY_H

#include "Utilities/PostProcessingStrategy/PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class FFTCalculator;
struct PhiAndNyTestParameterStruct;

class PhiAndNyTestPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
	static std::shared_ptr<PhiAndNyTestPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiAndNyTestParameterStruct> testPara);
	void evaluate();

	double getNyVx();
	double getNyVy();
	double getNyVz();
	double getPhiDiffVx();
	double getPhiDiffVy();
	double getPhiDiffVz();

private:
	PhiAndNyTestPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<PhiAndNyTestParameterStruct> testPara);
	
	std::vector<std::vector<double> > reduceDataToTimeSteps(std::vector<std::vector<double> > data, unsigned int startTimeStep, unsigned int endTimeStep);

	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::vector<std::string> dataToCalculatePhiAndNy;
	std::shared_ptr<FFTCalculator> fftCalculator;
	unsigned int startTimeStepCalculationPhiNy;
	unsigned int endTimeStepCalculationPhiNy;
	double nyVx, nyVy, nyVz;
	double phiDiffVx, phiDiffVy, phiDiffVz;
	bool isEvaluated;
};
#endif 