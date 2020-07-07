#ifndef L2_NORM_TEST_BETWEEN_KERNELS_H
#define L2_NORM_TEST_BETWEEN_KERNELS_H

#include "Utilities/Test/TestImp.h"

#include <memory>

class L2NormBetweenKernelPostProcessingStrategy;
class L2NormCalculator;
class L2NormCalculatorFactory;
class AnalyticalResults;
class SimulationResults;

class L2NormTestBetweenKernels : public TestImp
{
public:
	static std::shared_ptr<L2NormTestBetweenKernels> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep, std::string normalizeWith, std::shared_ptr<L2NormCalculatorFactory> factory);

	void update();
	void evaluate();
	std::string getLogFileOutput();
	std::string getErrorLogFileOutput();
	double getBasicL2Result();

	void setBasicSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy);
	void setDivergentKernelSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy);

private:
	L2NormTestBetweenKernels(std::shared_ptr<ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep, std::string normalizeWith, std::shared_ptr<L2NormCalculatorFactory> factory);
	int calcTimeStepInResults(unsigned int timeStep);
	std::vector<std::string> buildTestOutput();
	std::vector<std::string> buildBasicTestOutput();
	std::vector<std::string> buildErrorTestOutput();

	unsigned int timeStep;
	std::string dataToCalculate;
	std::shared_ptr<NumericalTestSimulation> basicSim;
	std::shared_ptr<SimulationInfo> basicSimInfo;
	std::shared_ptr<SimulationResults> basicSimResults;
	std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> basicPostProcessingStrategy;
	double basicL2Result;
	std::shared_ptr<NumericalTestSimulation> divergentSim;
	std::shared_ptr<SimulationInfo> divergentSimInfo;
	std::shared_ptr<SimulationResults> divergentSimResults;
	std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> divergentPostProcessingStrategy;
	double divergentL2Result;
	std::shared_ptr<L2NormCalculator> l2Normcalculator;
	std::string normalizeData;
	double resultL2ToBasicKernel;
};
#endif 