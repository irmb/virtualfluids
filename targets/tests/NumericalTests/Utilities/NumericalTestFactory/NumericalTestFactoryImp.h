#ifndef NUMERICAL_TEST_FACTORY_IMP_H
#define NUMERICAL_TEST_FACTORY_IMP_H

#include "NumericalTestFactory.h"
#include "Utilities\ConfigFileReader\ConfigData.h"

class AnalyticalResults;
class AnalyticalResults2DToVTKWriter;
class ColorConsoleOutput;
class L2NormTest;
class L2NormTestBetweenKernels;
class LogFileTimeInformation;
class LogFileQueueImp;
class PhiAndNuTest;
class PhiAndNuTestPostProcessingStrategy;
class SimulationInfo;
class SimulationLogFileInformation;
class SimulationParameter;
class SimulationResults;
class TestQueueImp;
class TestLogFileInformation;

class NumericalTestFactoryImp : public NumericalTestFactory
{
public:
	static std::shared_ptr< NumericalTestFactoryImp> getNewInstance(std::shared_ptr<ConfigDataStruct> configFileData);

	std::vector< std::shared_ptr< TestSimulation>> getTestSimulations();
	std::shared_ptr< TestQueue> getTestQueue();
	std::shared_ptr< LogFileQueue> getLogFileQueue();

private:
	NumericalTestFactoryImp() {};
	NumericalTestFactoryImp(std::shared_ptr<ConfigDataStruct> configFileData);

	void init();
	void makeTaylorGreenUxSimulations(std::string kernelName, double viscosity, double u0, double amplitude, int basicTimeStepLength);
	void makeTaylorGreenUzSimulations(std::string kernelName, double viscosity, double u0, double amplitude, int basicTimeStepLength);
	void makeShearWaveSimulations(std::string kernelName, double viscosity, double u0, double v0, int basicTimeStepLength);
	void makePeriodicBoundaryConditionSimulationAndTests(std::vector< std::shared_ptr< SimulationParameter>> simPara, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector< std::shared_ptr< AnalyticalResults>> analyResult, std::shared_ptr< SimulationLogFileInformation> simlogFileInfo, std::string kernelName, std::vector< bool> simulationsRun, double viscosity, int basicTimeStepLength);
	std::vector< std::shared_ptr< TestSimulation>> buildTestSimulation(std::vector< std::shared_ptr< SimulationParameter>> simPara, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector< std::shared_ptr< SimulationResults>> simResults, std::vector< std::shared_ptr< AnalyticalResults>> analyResult);
	std::vector< std::shared_ptr< PhiAndNuTest>> makePhiAndNuTests(std::vector<std::shared_ptr<TestSimulation>> testSim, std::vector<std::shared_ptr<SimulationInfo>> simInfo, std::vector< std::shared_ptr< PhiAndNuTestPostProcessingStrategy>> phiAndNuPostProStrategy, double viscosity, std::string dataToCalculate);
	std::vector< std::shared_ptr< L2NormTest>> makeL2NormTests(std::vector<std::shared_ptr< TestSimulation>> testSim, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector< std::shared_ptr< SimulationResults>> simResults, std::vector<std::shared_ptr<AnalyticalResults>> analyResult);
	std::vector< std::shared_ptr< L2NormTestBetweenKernels>> makeL2NormTestsBetweenKernels(std::vector<std::shared_ptr< TestSimulation>> testSim, std::vector< std::shared_ptr< SimulationInfo>> simInfo, std::vector< std::shared_ptr< SimulationResults>> simResults, std::vector<std::shared_ptr<AnalyticalResults>> analyResult);
	void makeLogFileWriter(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr< SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity, int basicTimeStepLength);
	
	void sortKernels();
	bool shouldSimulationGroupRun(std::vector<bool> test);
	unsigned int calcStartStepForToVectorWriter();
	bool checkNuAndPhiTestCouldRun(std::vector<bool> test);
	void calcNumberOfSimulations();

	std::shared_ptr<ConfigDataStruct> cfd;
	std::vector< std::shared_ptr< TestSimulation>> testSimulations;
	std::shared_ptr< TestQueueImp> testQueue;
	std::shared_ptr< LogFileQueueImp> logFileWriterQueue;
	std::shared_ptr< ColorConsoleOutput> colorOutput;
	std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;

	int simID;
	int numberOfSimulations;
	int simPerKernel, numberOfTestGroupsBetweenKernels, numberOfTestsForOneSimulation, numberOfTestsBetweenKernels;
	std::vector< std::shared_ptr< L2NormTestBetweenKernels>> l2KernelTests;
	int posBasicSimulationForL2Test, posDivergentSimulationForL2Test;
};
#endif