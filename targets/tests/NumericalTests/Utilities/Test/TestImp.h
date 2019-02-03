#ifndef TEST_IMP_H
#define TEST_IMP_H

#include "Utilities\Test\Test.h"

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
	void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PostProcessingStrategy> postProStrategy);

	virtual void evaluate() = 0;
	virtual std::vector< bool> getPassedTests() = 0;
	virtual void makeConsoleOutput() = 0;

	std::string getSimulationName();
	
protected:
	TestImp(std::shared_ptr< ColorConsoleOutput> colorOutput);
	bool CheckAllSimulationRun();

	std::vector<std::shared_ptr<NumericalTestSimulation> > simulations;
	std::vector<std::shared_ptr<PostProcessingStrategy> > postProStrategies;
	std::vector<std::shared_ptr<SimulationInfo> > simInfos;
	std::vector<bool> simulationRun;
	std::shared_ptr<ColorConsoleOutput> colorOutput;

	std::string kernelName;
	std::string simulationName;

private:
	TestImp() {};
};
#endif 