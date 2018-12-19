#ifndef TEST_SIMULATION_IMP_H
#define TEST_SIMULATION_IMP_H

#include "TestSimulation.h"

#include <vector>
#include <time.h>

class ToVectorWriter;
class ColorConsoleOutput;
class SimulationInfo;
class AnalyticalResults;
class AnalyticalResults2DToVTKWriter;

class TestSimulationImp : public TestSimulation
{
public:
	static std::shared_ptr< TestSimulation> getNewInsance(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< ColorConsoleOutput> colorOutput, std::shared_ptr< SimulationResults> simResults, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, bool writeAnalyticalResults);

	std::shared_ptr< SimulationParameter> getSimulationParameter();
	std::shared_ptr< DataWriter> getDataWriter();

	void registerSimulationObserver(std::shared_ptr< SimulationObserver> simObserver);
	bool getSimulationRun();
	std::shared_ptr< SimulationResults> getSimulationResults();

	void makeSimulationHeadOutput();
	void setSimulationStartTime();
	void setSimulationEndTimeAndNotifyObserver();
	std::string getRunTimeOutput();
	void setParameter(std::shared_ptr<Parameter> para);

private:
	TestSimulationImp(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< ColorConsoleOutput> colorOutput, std::shared_ptr< SimulationResults> simResults, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, bool writeAnalyticalResults);
	void notifyObserver();
	void setTestStartTime();
	void setTestEndTime();
	void setAnalyticalResultWriteTimeStartTime();
	void setAnalyticalResultWriteEndTime();
	double calcSimTime();
	float calcTestTime();
	double calcAnalyticalResultWriteTime();
	void writeAnalyticalResultsToVTK();

	std::shared_ptr< SimulationParameter> simPara;
	std::shared_ptr< SimulationInfo> simInfo;
	std::shared_ptr< SimulationResults> simResults;
	std::shared_ptr< ToVectorWriter> writeToVector;
	std::shared_ptr< ColorConsoleOutput> colorOutput;
	std::vector< std::shared_ptr< SimulationObserver>> simObserver;
	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;
	std::shared_ptr<Parameter> para;
	bool writeAnalyticalResults;

	bool simualtionRun;
	time_t simulationStartTime, simulationEndTime;
	clock_t testStartTime, testEndTime;
	time_t analyticalResultWriteStartTime, analyticalResultWriteEndTime;
	int simID;
};
#endif