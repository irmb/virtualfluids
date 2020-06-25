#ifndef TEST_SIMULATION_IMP_H
#define TEST_SIMULATION_IMP_H

#include "TestSimulation.h"
#include "Utilities/NumericalTestSimulation/NumericalTestSimulation.h"

#include <vector>
#include <time.h>

class ToVectorWriter;
class ColorConsoleOutput;
class SimulationInfo;
class AnalyticalResults;
class AnalyticalResults2DToVTKWriter;
class SimulationResults;
class TimeTracking;

struct TestSimulationDataStruct;

class TestSimulationImp : public TestSimulation, public NumericalTestSimulation
{
public:
	static std::shared_ptr<TestSimulationImp> getNewInsance(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResult, std::shared_ptr<TimeTracking> timeTracking, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput);

	std::shared_ptr<SimulationParameter> getSimulationParameter();
	std::shared_ptr<DataWriter> getDataWriter();
	std::shared_ptr<InitialCondition> getInitialCondition();
	std::shared_ptr<SimulationInfo> getSimulationInfo();
	std::shared_ptr<TimeTracking> getTimeTracking();

	SimulationStatus getSimulationStatus();
	
	void makeSimulationHeadOutput();
	void startPostProcessing();

	void setParameter(std::shared_ptr<Parameter> para);

	std::shared_ptr<SimulationResults> getSimulationResults();
	std::shared_ptr<AnalyticalResults> getAnalyticalResults();
	void registerSimulationObserver(std::shared_ptr<SimulationObserver> simObserver);
	std::vector<std::string> getDataToCalcTests();

private:
	TestSimulationImp(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResult, std::shared_ptr<TimeTracking> timeTracking, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput);
	void notifyObserver();

	void writeAnalyticalResultsToVTK();
	void checkSimulationResults();

	std::shared_ptr<SimulationParameter> simPara;
	std::shared_ptr<ToVectorWriter> toVectorWriter;
	std::shared_ptr<InitialCondition> initialCondition;
	std::shared_ptr<SimulationInfo> simInfo;
	std::shared_ptr<SimulationResults> simResult;
	std::shared_ptr<TimeTracking> timeTracking;

	std::shared_ptr<AnalyticalResults> analyticalResult;

	std::shared_ptr<ColorConsoleOutput> colorOutput;
	std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;
	std::shared_ptr<Parameter> para;
	std::vector<std::shared_ptr<SimulationObserver> > simObserver;
	
	std::vector<std::string> dataToCalcTests;
	SimulationStatus status;
};
#endif