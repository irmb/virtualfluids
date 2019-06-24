#ifndef LOGFILE_DATA_ASSISTANT_IMP_H
#define LOGFILE_DATA_ASSISTANT_IMP_H

#include "LogFileDataAssistant.h"

class LogFileDataGroupImp;
class LogFileDataAssistantStrategy;

class LogFileDataAssistantImp : public LogFileDataAssistant
{
public:
	static std::shared_ptr<LogFileDataAssistant> getNewInstance();


	std::vector<std::shared_ptr<LogFileDataGroup> > findDataCombination(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy, DataCombination combination);
	

protected:
	LogFileDataAssistantImp();
	
	std::vector<std::shared_ptr<LogFileDataGroup> > findEqualSimulationsForDifferentKernels(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy);
	std::vector<std::shared_ptr<LogFileDataGroup> > findEqualKernelSimulationsForDifferentViscosities(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy);

	std::vector<std::shared_ptr<LogFileData> > getSimulationGroupLogFileData(std::string simName, std::vector<std::shared_ptr<LogFileData> > allLogFileData);
	std::vector<std::vector<std::shared_ptr<LogFileData> > > sortLogFileDataAfterKernels(std::vector<std::shared_ptr<LogFileData> > logFileData);

	bool checkEqualSimulationData(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);
	bool checkEqualViscosity(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);
	bool checkEqualKernel(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);

	bool checkBasicSimulationIsInLogFiles(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::string simName);

	std::vector<std::shared_ptr<LogFileDataGroup>> castToLogFileDataGroup(std::vector<std::shared_ptr<LogFileDataGroupImp>> data);

	bool equalDouble(double num1, double num2);

};
#endif