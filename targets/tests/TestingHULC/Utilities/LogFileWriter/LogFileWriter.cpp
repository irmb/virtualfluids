#include "LogFileWriter.h"

#include <helper_functions.h>

#include <iomanip>
#include <ctime>

LogFileWriter::LogFileWriter(std::string filePath)
{
	std::ostringstream oss;
	oss << filePath << "\\logFile_" << calcDateAndTime() << ".txt";
	this->logFilePath = oss.str();

	logFile.open(logFilePath, std::ios::out);
}

std::shared_ptr<LogFileWriter> LogFileWriter::getNewInstance(std::string filePath)
{
	return std::shared_ptr<LogFileWriter>(new LogFileWriter(filePath));
}

void LogFileWriter::makeOutput(std::string output)
{
	logFile << output;
}

std::string LogFileWriter::calcDateAndTime()
{
	std::ostringstream oss;
	now = time(NULL);
	nowLocal = *localtime(&now);
	oss << std::setfill('0')  << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
	return oss.str();
}