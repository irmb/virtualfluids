#ifndef L2_NORM_TEST_H
#define L2_NORM_TEST_H

#include "Utilities\Test\TestImp.h"

#include <memory>

class L2NormCalculator;
class AnalyticalResults;

class L2NormTest : public TestImp
{
public:
	static std::shared_ptr<L2NormTest> getNewInstance(std::shared_ptr< AnalyticalResults> analyticalResult, std::shared_ptr< ColorConsoleOutput> colorOutput);

	void update();
	void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo);
	void evaluate();
	std::string getLogFileOutput();
	std::vector< bool> getPassedTests();
	void makeConsoleOutput();

private:
	L2NormTest(std::shared_ptr< AnalyticalResults> analyticalResult, std::shared_ptr< ColorConsoleOutput> colorOutput);

	std::shared_ptr< L2NormCalculator> calculator;
	std::shared_ptr< AnalyticalResults> analyticalResult;
};
#endif 