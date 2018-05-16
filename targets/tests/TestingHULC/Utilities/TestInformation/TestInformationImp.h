#ifndef TEST_INFORMATION_IMP_H
#define TEST_INFORMATION_IMP_H

#include "TestInformation.h"

#include <sstream>
#include <memory>
#include <vector>

class SimulationInfo;
class LogFileInformation;

class TestInformationImp : public TestInformation
{
public:
	void makeSimulationHeadOutput(int i);
	void setSimulationStartTime(int i);
	void setSimulationEndTime(int i);
	void writeLogFile();

	static std::shared_ptr<TestInformationImp> getNewInstance();
	void setSimulationInfo(std::vector< std::shared_ptr<SimulationInfo> > simInfos);
	void setLogFilePath(std::string aLogFilePath);
	void setLogFileInformation(std::vector< std::shared_ptr<LogFileInformation> > logInfos);

protected:
	TestInformationImp();

private:
	void makeHashLine();
	void makeCenterHead(std::string head);

	std::vector <std::shared_ptr<SimulationInfo> > simInfos;
	std::vector< std::shared_ptr<LogFileInformation> > logInfos;
	int numberOfTimeSteps;
	int basisTimeStepLength;
	int startStepCalculation;
	double viscosity;
	
	std::ostringstream oss;
	std::string logFilePath;

};
#endif 
