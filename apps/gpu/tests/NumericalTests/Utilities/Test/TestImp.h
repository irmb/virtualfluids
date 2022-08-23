#ifndef TEST_IMP_H
#define TEST_IMP_H

#include "Utilities/Test/Test.h"

#include <vector>
#include <memory>

class NumericalTestSimulation;
class SimulationResults;
class SimulationInfo;
class ColorConsoleOutput;
class PostProcessingStrategy;

class TestImp : public Test
{
public:
	void run() override;
	void update();
	TestStatus getTestStatus();
	virtual void makeConsoleOutput();

	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PostProcessingStrategy> postProStrategy);
		
protected:
	TestImp(std::shared_ptr<ColorConsoleOutput> colorOutput);

	virtual void evaluate() = 0;
	virtual std::vector<std::string> buildTestOutput() = 0;
	virtual std::vector<std::string> buildBasicTestOutput() = 0;
	virtual std::vector<std::string> buildErrorTestOutput() = 0;
	std::vector<std::string> buildSimulationFailedTestOutput();
	
	bool CheckAllSimulationRun();

	std::vector<std::shared_ptr<NumericalTestSimulation> > simulations;
	std::vector<std::shared_ptr<PostProcessingStrategy> > postProStrategies;
	std::vector<std::shared_ptr<SimulationInfo> > simInfos;
	std::vector<bool> simulationRun;
	std::shared_ptr<ColorConsoleOutput> colorOutput;
	TestStatus testStatus;

private:
	TestImp() {};
};
#endif 