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
	void update();
	void setSimulationCrashed();
	TestStatus getTestStatus();
	virtual void makeConsoleOutput();

	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PostProcessingStrategy> postProStrategy);
		
protected:
	virtual void evaluate() = 0;
	virtual std::vector<std::string> buildTestOutput() = 0;
	virtual std::vector<std::string> buildBasicTestOutput() = 0;
	virtual std::vector<std::string> buildErrorTestOutput() = 0;

	TestImp(std::shared_ptr<ColorConsoleOutput> colorOutput);
	bool CheckAllSimulationRun();


	std::vector<std::string> buildSimulationFailedTestOutput();
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