#ifndef LOG_FILE_WRITER_IMP_H
#define LOG_FILE_WRITER_IMP_H

#include "LogFileWriter.h"

#include <fstream>
#include <vector>
#include <memory>

class SimulationLogFileInformation;
class LogFileInformation;
class LogFileTimeInformation;
class TestLogFileInformation;

class LogFileWriterImp : public LogFileWriter
{
public:
	static std::shared_ptr<LogFileWriterImp> getNewInstance(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles,
														std::shared_ptr< LogFileTimeInformation> logFileTimeInfo,
														std::shared_ptr< SimulationLogFileInformation> simLogInfo, 
														std::string kernelName, 
														double viscosity, std::vector<int> devices, int numberOfTimeSteps, int startStepCalculation, int basicTimeStepLength);
	void writeLogFile(std::string basicFilePath);
	

private:
	LogFileWriterImp(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr< SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity, std::vector<int> devices, int numberOfTimeSteps, int startStepCalculation, int basicTimeStepLength);
	std::string calcDateAndTime();
	std::string buildFilePath(std::string basicFilePath);

	std::fstream logFile;
	std::string logFilePath;
	time_t now;
	struct tm nowLocal;
	std::string kernelName;
	double viscosity;
	std::vector< std::shared_ptr< LogFileInformation>> logFileInfo;
	std::shared_ptr< SimulationLogFileInformation> simLogInfo;
};
#endif 
