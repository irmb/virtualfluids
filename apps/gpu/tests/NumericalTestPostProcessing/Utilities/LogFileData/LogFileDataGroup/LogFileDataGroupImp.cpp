#include "LogFileDataGroupImp.h"

std::shared_ptr<LogFileDataGroupImp> LogFileDataGroupImp::getNewInstance()
{
	return std::shared_ptr<LogFileDataGroupImp>(new LogFileDataGroupImp());
}

std::shared_ptr<LogFileData> LogFileDataGroupImp::getLogFileData(int number)
{
	return data.at(number);
}

int LogFileDataGroupImp::getGroupSize()
{
	return data.size();
}

void LogFileDataGroupImp::addLogFileData(std::shared_ptr<LogFileData> logFileData)
{
	data.push_back(logFileData);
}

LogFileDataGroupImp::LogFileDataGroupImp()
{
	data.resize(0);
}
