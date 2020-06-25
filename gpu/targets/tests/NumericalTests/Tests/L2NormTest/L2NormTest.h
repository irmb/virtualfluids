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
	static std::shared_ptr<L2NormTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData);

	void update();
	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<L2NormPostProcessingStrategy> postProStrategy);

	void evaluate();
	std::string getLogFileOutput();
	std::string getErrorLogFileOutput();

private:
	L2NormTest(std::shared_ptr<ColorConsoleOutput> colorOutput, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::string dataToCalculate, double maxL2NormDiff, std::string normalizeData);
	std::vector<std::string> buildTestOutput();
	std::vector<std::string> buildBasicTestOutput();
	std::vector<std::string> buildErrorTestOutput();


	unsigned int basicTimeStep, divergentTimeStep;
	double resultBasicTimestep, resultDivergentTimeStep;
	std::string dataToCalculate;
	double diffL2Norm;
	double maxL2NormDiff;
	bool testPassed;

	std::string normalizeData;

	std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > l2NormPostProStrategies;
};
#endif 