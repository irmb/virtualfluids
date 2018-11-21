#ifndef CONFIG_FILE_READER_H
#define CONFIG_FILE_READER_H

#include "LBM\LB.h"

#include <memory>
#include <vector>
#include <string>

class SimulationParameter;
class SimulationInfo;
class TestCout;
class PhiAndNuTest;
class TestSimulation;
class TestQueueImp;
class TestQueue;
class LogFileQueueImp;
class LogFileQueue;
class LogFileInformation;
class LogFileTimeInformation;
class TestLogFileInformation;
class SimulationLogFileInformation;

class ConfigFileReader
{
public:
	static std::shared_ptr< ConfigFileReader> getNewInstance();

	void readConfigFile(const std::string aFilePath);

	std::vector< std::shared_ptr< TestSimulation>> getTestSimulations();
	std::shared_ptr< TestQueue> getTestQueue();
	std::shared_ptr< LogFileQueue> getLogFileQueue();

protected:
	ConfigFileReader();
	
private:
	void checkConfigFileData();
	void init();

	void makeLogFileWriter(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr< SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity);

	std::vector< std::shared_ptr< TestSimulation>> buildTestSimulation(std::vector< std::shared_ptr< SimulationParameter>> simPara, std::vector< std::shared_ptr< SimulationInfo>> simInfo);
	void makeTaylorGreenSimulations(std::string kernelName, double viscosity, double u0, double amplitude);
	void makeShearWaveSimulations(std::string kernelName, double viscosity, double u0, double v0);

	std::vector< std::shared_ptr< PhiAndNuTest>> makePhiAndNuTests(std::vector< std::shared_ptr< TestSimulation>> testSim, std::vector< std::shared_ptr< SimulationInfo>> simInfo, double viscosity);
	
	bool shouldSimulationGroupRun(std::vector<bool> test);


	std::vector<double> u0SW, v0SW;
	std::vector<double> amplitudeTGV, u0TGV;
	bool nuAndPhiTestTGV, nuAndPhiTestSW;
	std::string dataToCalcPhiAndNuTest;
	std::vector<double> viscosity;
	real rho0;
	real l0;
	double minOrderOfAccuracy;
	unsigned int numberOfTimeSteps, basisTimeStepLength, startStepCalculation;
	unsigned int ySliceForCalculation;
	unsigned int startStepFileWriter;
	unsigned int maxLevel;
	unsigned int numberOfGridLevels;
	bool writeFiles;
	std::string filePath;
	std::string logFilePath;
	std::vector< std::string> kernelsToTest;
	std::vector< std::string> grids;
	std::vector< real> lx;
	std::vector< real> lz;
	std::vector< int> devices;
	std::vector< bool> tgv;
	std::vector< bool> sw;
	
	std::vector< std::shared_ptr< LogFileInformation> > logInfo;

	std::shared_ptr< TestQueueImp> testQueue;
	std::vector< std::shared_ptr< TestSimulation>> testSimulation;

	std::shared_ptr< LogFileQueueImp> logFileWriterQueue;

	int simID;
};
#endif