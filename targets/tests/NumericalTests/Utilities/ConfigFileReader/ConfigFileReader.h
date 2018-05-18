#ifndef CONFIG_FILE_READER_H
#define CONFIG_FILE_READER_H

#include "LBM\LB.h"

#include <memory>
#include <vector>
#include <string>

class TestInformation;
class TestParameter;
class TestResults;
class TestCout;
class TestInformationImp;
class PhiAndNuTest;
class SimulationInfo;
class LogFileInformation;

class ConfigFileReader
{
public:
	static std::shared_ptr< ConfigFileReader> getNewInstance();

	void readConfigFile(const std::string aFilePath);
	std::shared_ptr< TestInformation> getTestInformation();
	std::vector< std::shared_ptr< TestParameter> > getTestParameter();

protected:
	ConfigFileReader();
	
private:
	void makeTestInformation();
	void makeTestParameter();
	void makeSimulationInfo();
	void makeLogFileInformation();
	void makeTestResults();
	bool testShouldRun(std::vector<bool> test);

	real viscosity, rho0;
	real u0SW, v0SW;
	real amplitudeTGV, u0TGV;
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
	std::vector< std::string> grids;
	std::vector< real> lx;
	std::vector< real> lz;
	std::vector< int> devices;
	std::vector< bool> tgv;
	std::vector< bool> sw;
	


	std::shared_ptr< TestCout> testOutput;
	std::shared_ptr< TestInformationImp> testInfo;
	std::vector< std::shared_ptr< PhiAndNuTest> > tests;
	std::vector< std::shared_ptr< TestParameter> > testParameter;
	std::vector< std::shared_ptr< SimulationInfo> > simInfo;
	std::vector< std::shared_ptr< LogFileInformation> > logInfo;
	std::vector< std::shared_ptr< TestResults> > testResults;
};
#endif