#ifndef TEST_SIMULATION_IMP_H
#define TEST_SIMULATION_IMP_H

#include "TestSimulation.h"

#include <vector>
#include <time.h>

class ToVectorWriter;
class ColorConsoleOutput;
class SimulationInfo;

class TestSimulationImp : public TestSimulation
{
public:
	static std::shared_ptr< TestSimulation> getNewInsance(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< ColorConsoleOutput> colorOutput);

	std::shared_ptr< SimulationParameter> getSimulationParameter();
	std::shared_ptr< DataWriter> getDataWriter();

	void registerSimulationObserver(std::shared_ptr< SimulationObserver> simObserver);
	bool getSimulationRun();
	std::shared_ptr< SimulationResults> getSimulationResults();

	void makeSimulationHeadOutput();
	void setSimulationStartTime();
	void setSimulationEndTimeAndNotifyObserver();
	std::string getRunTimeOutput();

private:
	TestSimulationImp(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< ColorConsoleOutput> colorOutput);
	void notifyObserver();
	void setTestStartTime();
	void setTestEndTime();
	double calcSimTime();
	float calcTestTime();

	std::shared_ptr< SimulationParameter> simPara;
	std::shared_ptr< SimulationInfo> simInfo;
	std::shared_ptr< SimulationResults> simResults;
	std::shared_ptr< ToVectorWriter> writeToVector;
	std::shared_ptr< ColorConsoleOutput> colorOutput;
	std::vector< std::shared_ptr< SimulationObserver>> simObserver;

	bool simualtionRun;
	time_t simulationStartTime, simulationEndTime;
	clock_t testStartTime, testEndTime;
	int simID;
};
#endif