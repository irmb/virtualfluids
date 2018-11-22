#include "LogFileWriterImp.h"

#include "Utilities\LogFileInformation\SimulationLogFileInformation\SimulationLogFileInformation.h"
#include "Utilities\LogFileInformation\LogFileHead\LogFileHead.h"
#include "Utilities\LogFileInformation\BasicSimulationInfo\BasicSimulationInfo.h"
#include "Utilities\LogFileInformation\LogFileTimeInformation\LogFileTimeInformation.h"
#include "Utilities\LogFileInformation\TestLogFileInformation\TestLogFileInformation.h"

#include <helper_functions.h>
#include <iomanip>
#include <ctime>
#include <experimental/filesystem>

LogFileWriterImp::LogFileWriterImp(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr< SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity, std::vector<int> devices, int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation) : kernelName(kernelName), viscosity(viscosity)
{
	logFileInfo.push_back(LogFileHead::getNewInstance(devices));
	logFileInfo.push_back(BasicSimulationInfo::getNewInstance(numberOfTimeSteps, basisTimeStepLength, viscosity));
	this->simLogInfo = simLogInfo;
	logFileInfo.push_back(this->simLogInfo);
	logFileInfo.push_back(logFileTimeInfo);
	for (int i = 0; i < testLogFiles.size(); i++)
		logFileInfo.push_back(testLogFiles.at(i));
}

std::shared_ptr<LogFileWriterImp> LogFileWriterImp::getNewInstance(std::vector< std::shared_ptr< TestLogFileInformation>> testLogFiles, std::shared_ptr< LogFileTimeInformation> logFileTimeInfo, std::shared_ptr< SimulationLogFileInformation> simLogInfo, std::string kernelName, double viscosity, std::vector<int> devices, int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation)
{
	return std::shared_ptr<LogFileWriterImp>(new LogFileWriterImp(testLogFiles, logFileTimeInfo, simLogInfo, kernelName, viscosity, devices, numberOfTimeSteps, basisTimeStepLength, startStepCalculation));
}

void LogFileWriterImp::writeLogFile(std::string basicFilePath)
{
	logFilePath = buildFilePath(basicFilePath);

	logFile.open(logFilePath, std::ios::out);

	bool test = logFile.is_open();

	for (int i = 0; i < logFileInfo.size(); i++)
		logFile << logFileInfo.at(i)->getOutput();	
}


std::string LogFileWriterImp::calcDateAndTime()
{
	std::ostringstream oss;
	now = time(NULL);
	nowLocal = *localtime(&now);
	oss << std::setfill('0')  << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
	return oss.str();
}

std::string LogFileWriterImp::buildFilePath(std::string basicFilePath)
{
	std::ostringstream filePath;
	filePath << basicFilePath << kernelName << "\\viscosity_" << viscosity << "\\" << simLogInfo->getFilePathExtension();
	
	std::experimental::filesystem::path dir(filePath.str());
	if (!(std::experimental::filesystem::exists(dir)))
		std::experimental::filesystem::create_directories(dir);

	filePath << "\\logfile_" << calcDateAndTime() << "_" << kernelName << "_vis_" << viscosity << ".txt";
	return filePath.str();
}
