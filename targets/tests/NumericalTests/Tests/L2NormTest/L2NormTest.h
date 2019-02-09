#ifndef L2_NORM_TEST_H
#define L2_NORM_TEST_H

#include "Utilities/Test/TestImp.h"

#include <memory>

class L2NormCalculator;
class AnalyticalResults;
class L2NormPostProcessingStrategy;

struct L2NormTestParameterStruct;

class L2NormTest : public TestImp
{
public:
	static std::shared_ptr<L2NormTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate);

	void update();
	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormPostProcessingStrategy> postProStrategy);
	void evaluate();
	std::string getLogFileOutput();
	std::vector<bool> getPassedTests();
	void makeConsoleOutput();

private:
	L2NormTest(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate);

	unsigned int basicTimeStep, divergentTimeStep;
	double resultBasicTimestep, resultDivergentTimeStep;
	std::string dataToCalculate;
	double diffL2Norm;
	double maxL2NormDiff;
	bool testPassed;

	std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > l2NormPostProStrategies;
};
#endif 