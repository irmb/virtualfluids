#ifndef LOGFILE_DATA_ASSISTANT_IMP_H
#define LOGFILE_DATA_ASSISTANT_IMP_H

#include "LogFileDataAssistant.h"

class LogFileDataAssistantStrategyFactory;
class LogFileDataGroupImp;

class LogFileDataAssistantImp : public LogFileDataAssistant
{
public:
	static std::shared_ptr<LogFileDataAssistant> getNewInstance();


	std::vector<std::shared_ptr<LogFileDataGroup> > findDataCombination(std::vector<std::shared_ptr<LogFileData> > allLogFileData, BasicSimulation simulation, DataCombination combination);
	


protected:
	LogFileDataAssistantImp();
	
	std::vector<std::shared_ptr<LogFileDataGroup> > findEqualSimulationsForDifferentKernels(std::vector<std::shared_ptr<LogFileData> > allLogFileData, BasicSimulation simulation);
	std::vector<std::shared_ptr<LogFileDataGroup> > findEqualKernelSimulationsForDifferentViscosities(std::vector<std::shared_ptr<LogFileData> > allLogFileData, BasicSimulation simulation);

	std::vector<std::shared_ptr<LogFileData> > getSimulationGroupLogFileData(std::string simName, std::vector<std::shared_ptr<LogFileData> > allLogFileData);
	std::vector<std::vector<std::shared_ptr<LogFileData> > > sortLogFileDataAfterKernels(std::vector<std::shared_ptr<LogFileData> > logFileData);

	bool checkEqualSimulationData(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);
	bool checkEqualViscosity(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);
	bool checkEqualKernel(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);

	std::vector<std::shared_ptr<LogFileDataGroup>> castToLogFileDataGroup(std::vector<std::shared_ptr<LogFileDataGroupImp>> data);

	bool equalDouble(double num1, double num2);

	std::shared_ptr<LogFileDataAssistantStrategyFactory> assistentStrategyFactory;
};
#endif