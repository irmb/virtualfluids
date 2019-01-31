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

struct TestSimulationDataStruct;

class TestSimulationImp : public TestSimulation
{
public:
	static std::shared_ptr< TestSimulationImp> getNewInsance(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResults, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput);

	std::shared_ptr< SimulationParameter> getSimulationParameter();
	std::shared_ptr< DataWriter> getDataWriter();
	std::shared_ptr<SimulationInfo> getSimulationInfo();

	void registerSimulationObserver(std::shared_ptr< SimulationObserver> simObserver);
	bool getSimulationRun();
	std::shared_ptr< SimulationResults> getSimulationResults();
	std::shared_ptr<AnalyticalResults> getAnalyticalResults();

	void makeSimulationHeadOutput();
	void setSimulationStartTime();
	void setSimulationEndTimeAndNotifyObserver();
	std::string getRunTimeOutput();
	void setParameter(std::shared_ptr<Parameter> para);

private:
	TestSimulationImp(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResults, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput);
	void notifyObserver();
	void setTestStartTime();
	void setTestEndTime();
	void setAnalyticalResultWriteTimeStartTime();
	void setAnalyticalResultWriteEndTime();
	double calcSimTime();
	float calcTestTime();
	double calcAnalyticalResultWriteTime();
	void writeAnalyticalResultsToVTK();

	std::shared_ptr<SimulationParameter> simPara;
	std::shared_ptr<SimulationInfo> simInfo;
	std::shared_ptr<SimulationResults> simResults;
	std::shared_ptr<ToVectorWriter> toVectorWriter;
	std::shared_ptr<ColorConsoleOutput> colorOutput;
	std::shared_ptr<AnalyticalResults> analyticalResult;
	std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;
	std::shared_ptr<Parameter> para;
	std::vector< std::shared_ptr< SimulationObserver>> simObserver;
	bool simualtionRun;
	time_t simulationStartTime, simulationEndTime;
	clock_t testStartTime, testEndTime;
	time_t analyticalResultWriteStartTime, analyticalResultWriteEndTime;
};
#endif