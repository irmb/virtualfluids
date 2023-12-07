#include "LogFileQueueImp.h"

#include <helper_functions.h>

#include <ctime>
#include <iomanip>

#include <basics/DataTypes.h>

#include "Utilities/LogFileWriter/LogFileWriter.h"

std::shared_ptr<LogFileQueueImp> LogFileQueueImp::getNewInstance(std::string basicLogFilePath)
{
    return std::shared_ptr<LogFileQueueImp>(new LogFileQueueImp(basicLogFilePath));
}

void LogFileQueueImp::writeLogFiles()
{
    for (uint i = 0; i < logFileWriter.size(); i++){
        logFileWriter.at(i)->writeLogFile(basicLogFilePath);
    }
}

void LogFileQueueImp::addLogFileWriter(std::shared_ptr<LogFileWriter> aLogFileWriter)
{
    logFileWriter.push_back(aLogFileWriter);
}

LogFileQueueImp::LogFileQueueImp(std::string basicLogFilePath)
{
    logFileWriter.resize(0);

    std::ostringstream oss;
    oss << basicLogFilePath << "/NumericalTestLogFiles/";
    this->basicLogFilePath = oss.str();
}

std::string LogFileQueueImp::calcDateAndTime()
{
    std::ostringstream oss;
    now = time(NULL);
    nowLocal = *localtime(&now);
    oss << std::setfill('0') << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
    return oss.str();
}
