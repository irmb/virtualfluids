#ifndef L2_NORM_TEST_BETWEEN_KERNELS_H
#define L2_NORM_TEST_BETWEEN_KERNELS_H

#include "Utilities\Test\TestImp.h"

#include <memory>

class L2NormBetweenKernelPostProcessingStrategy;
class L2NormCalculator;
class AnalyticalResults;
class SimulationResults;

class L2NormTestBetweenKernels : public TestImp
{
public:
	static std::shared_ptr<L2NormTestBetweenKernels> getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep);

	void update();
	void evaluate();
	std::string getLogFileOutput();
	double getBasicL2Result();
	std::vector< bool> getPassedTests();
	void makeConsoleOutput();

	void setBasicSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy);
	void setDivergentKernelSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< L2NormBetweenKernelPostProcessingStrategy> postProcessingStrategy);

private:
	L2NormTestBetweenKernels(std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, unsigned int timeStep);
	int calcTimeStepInResults(unsigned int timeStep);

	unsigned int timeStep;
	std::string dataToCalculate;
	bool testPassed;

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

	double resultL2ToBasicKernel;
};
#endif 