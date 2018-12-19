#ifndef PHIANDNU_TEST_POST_PROCESSING_STRATEGY_H
#define PHIANDNU_TEST_POST_PROCESSING_STRATEGY_H

#include "Utilities\PostProcessingStrategy\PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class FFTCalculator;

class PhiAndNuTestPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
	static std::shared_ptr< PhiAndNuTestPostProcessingStrategy> getNewInstance(std::shared_ptr< SimulationResults> simResult, std::shared_ptr< AnalyticalResults> analyticalResult, std::vector<std::string> dataToCalculatePhiAndNu , unsigned int startTimeStepCalculationPhiNu, unsigned int endTimeStepCalculationPhiNu);
	void evaluate();

	double getNuVx();
	double getNuVy();
	double getNuVz();
	double getPhiDiffVx();
	double getPhiDiffVy();
	double getPhiDiffVz();

private:
	PhiAndNuTestPostProcessingStrategy(std::shared_ptr< SimulationResults> simResult, std::shared_ptr< AnalyticalResults> analyticalResult, std::vector<std::string> dataToCalculatePhiAndNu, unsigned int startTimeStepCalculationPhiNu, unsigned int endTimeStepCalculationPhiNu);
	
	std::vector<std::vector<double>> reduceDataToTimeSteps(std::vector<std::vector<double>> data, unsigned int startTimeStep, unsigned int endTimeStep);

	std::shared_ptr< AnalyticalResults> analyticalResult;
	std::vector<std::string> dataToCalculatePhiAndNu;
	std::shared_ptr< FFTCalculator> fftCalculator;
	unsigned int startTimeStepCalculationPhiNu;
	unsigned int endTimeStepCalculationPhiNu;
	double nuVx, nuVy, nuVz;
	double phiDiffVx, phiDiffVy, phiDiffVz;
	bool isEvaluated;
};
#endif 