#ifndef L2NORM_BETWEEN_KERNEL_POST_PROCESSING_STRATEGY_H
#define L2NORM_BETWEEN_KERNEL_POST_PROCESSING_STRATEGY_H

#include "Utilities\PostProcessingStrategy\PostProcessingStrategyImp.h"

#include <memory>

class AnalyticalResults;
class L2NormCalculator;
struct L2NormTestBetweenKernelsParameterStruct;

class L2NormBetweenKernelPostProcessingStrategy : public PostProcessingStrategyImp
{
public:
	static std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> getNewInstance(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara);
	void evaluate();

	double getL2NormVx(int timeStep);
	double getL2NormVy(int timeStep);
	double getL2NormVz(int timeStep);
	double getL2NormPress(int timeStep);
	double getL2NormRho(int timeStep);

	virtual std::shared_ptr<SimulationResults> getSimulationResult();

private:
	L2NormBetweenKernelPostProcessingStrategy(std::shared_ptr<SimulationResults> simResult, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara);

	int calcPosInTimeStep(int time);

	std::vector<int> timeSteps;
	std::vector<double> l2Norm;
	std::vector<double> l2Vx;
	std::vector<double> l2Vy;
	std::vector<double> l2Vz;
	std::vector<double> l2Press;
	std::vector<double> l2Rho;
	std::vector<std::string> dataToCalculate;
	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::shared_ptr<L2NormCalculator> l2Normcalculator;
	bool isEvaluated;

};
#endif