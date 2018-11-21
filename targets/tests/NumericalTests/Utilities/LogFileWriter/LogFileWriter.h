#ifndef LOG_FILE_WRITER_H
#define LOG_FILE_WRITER_H

#include <fstream>
#include <vector>
#include <memory>
#include <string>


class SimulationLogFileInformation;
class LogFileInformation;
class LogFileTimeInformation;
class TestLogFileInformation;

class LogFileWriter
{
public:
	static std::shared_ptr<LogFileWriter> getNewInstance(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, 
														std::shared_ptr< LogFileTimeInformation> logFileTimeInfo,
														std::shared_ptr< SimulationLogFileInformation> simLogInfo, 
														std::string kernelName, 
														double viscosity, std::vector<int> devices, int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation);
	void writeLogFile(std::string basicFilePath);
	

private:
	LogFileWriter(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr< SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity, std::vector<int> devices, int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation);
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
