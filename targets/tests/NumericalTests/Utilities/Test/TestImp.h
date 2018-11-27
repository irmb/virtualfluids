#ifndef TEST_IMP_H
#define TEST_IMP_H

#include "Utilities\Test\Test.h"

#include <vector>

class TestSimulation;
class SimulationResults;
class SimulationInfo;
class ColorConsoleOutput;

class TestImp : public Test
{
public:
	void update();
	void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< PostProcessingResults> postProResult);

	virtual void evaluate() = 0;
	virtual std::string getLogFileOutput() = 0;
	virtual std::vector< bool> getPassedTests() = 0;
	virtual void makeConsoleOutput() = 0;

	std::string getSimulationName();
	
protected:
	TestImp(std::shared_ptr< ColorConsoleOutput> colorOutput);
	bool CheckAllSimulationRun();


	std::vector< bool> simulationRun;
	std::shared_ptr< ColorConsoleOutput> colorOutput;
	std::vector< std::shared_ptr< TestSimulation>> simulations;
	std::vector< std::shared_ptr< PostProcessingResults>> postProResults;
	std::vector< std::shared_ptr< SimulationInfo>> simInfos;

	std::string kernelName;
	std::string simulationName;

private:
	TestImp() {};
};
#endif 