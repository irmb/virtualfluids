#ifndef L2_NORM_TEST_H
#define L2_NORM_TEST_H

#include "Utilities\Test\TestImp.h"

#include <memory>

class L2NormCalculator;
class AnalyticalResults;
class PostProcessingResults;

class L2NormTest : public TestImp
{
public:
	static std::shared_ptr<L2NormTest> getNewInstance(std::shared_ptr< PostProcessingResults> postProResults, std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double maxL2NormDiff, unsigned int basicTimeStep, unsigned int divergentTimeStep);

	void update();
	void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< PostProcessingResults> postProResults);
	void evaluate();
	std::string getLogFileOutput();
	std::vector< bool> getPassedTests();
	void makeConsoleOutput();

private:
	L2NormTest(std::shared_ptr< PostProcessingResults> postProResults, std::shared_ptr< ColorConsoleOutput> colorOutput, std::string dataToCalculate, double maxL2NormDiff, unsigned int basicTimeStep, unsigned int divergentTimeStep);

	std::shared_ptr< PostProcessingResults> postProResults;
	unsigned int basicTimeStep, divergentTimeStep;
	double resultBasicTimestep, resultDivergentTimeStep;
	std::string dataToCalculate;
	double diffL2Norm;
	double maxL2NormDiff;
	bool testPassed;
};
#endif 