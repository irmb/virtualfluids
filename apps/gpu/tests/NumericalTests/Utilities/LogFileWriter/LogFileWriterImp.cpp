#include "LogFileWriterImp.h"

#include "Utilities/LogFileInformation/SimulationLogFileInformation/SimulationLogFileInformation.h"
#include "Utilities/LogFileInformation/LogFileHead/LogFileHead.h"
#include "Utilities/LogFileInformation/BasicSimulationInfo/BasicSimulationInfo.h"
#include "Utilities/LogFileInformation/BasicTestLogFileInformation/BasicTestLogFileInformation.h"
#include "Utilities/LogFileInformation/LogFileTimeInformation/LogFileTimeInformation.h"
#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"

#include <helper_functions.h>
#include <iomanip>
#include <ctime>
#include <filesystem>

LogFileWriterImp::LogFileWriterImp(std::shared_ptr<LogFileHead> logFileHead, std::shared_ptr<BasicSimulationInfo> basicSimInfo, std::shared_ptr<BasicTestLogFileInformation> basicTestInfo, std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<LogFileTimeInformation> logFileTimeInfo, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::string kernel, double viscosity) : viscosity(viscosity)
{
    kernelName = kernel;

    logFileInfo.push_back(logFileHead);
    logFileInfo.push_back(basicSimInfo);
    this->simLogInfo = simLogInfo;
    logFileInfo.push_back(simLogInfo);
    logFileInfo.push_back(logFileTimeInfo);
    logFileInfo.push_back(basicTestInfo);
    for (int i = 0; i < testLogFiles.size(); i++)
        logFileInfo.push_back(testLogFiles.at(i));
}

std::shared_ptr<LogFileWriterImp> LogFileWriterImp::getNewInstance(std::shared_ptr<LogFileHead> logFileHead, std::shared_ptr<BasicSimulationInfo> basicSimInfo, std::shared_ptr<BasicTestLogFileInformation> basicTestInfo, std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<LogFileTimeInformation> logFileTimeInfo, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::string kernel, double viscosity)
{
    return std::shared_ptr<LogFileWriterImp>(new LogFileWriterImp(logFileHead, basicSimInfo, basicTestInfo, testLogFiles, logFileTimeInfo, simLogInfo, kernel, viscosity));
}

void LogFileWriterImp::writeLogFile(std::string basicFilePath)
{
    logFilePath = buildFilePath(basicFilePath);
    logFile.open(logFilePath, std::ios::out);

    for (int i = 0; i < logFileInfo.size(); i++)
        logFile << logFileInfo.at(i)->getOutput();    

    logFile.close();
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
    filePath << basicFilePath << simLogInfo->getFilePathExtension().at(0) << "/viscosity_" << viscosity << "/" << simLogInfo->getFilePathExtension().at(1) << "/" << kernelName;
    
    std::filesystem::path dir(filePath.str());
    if (!(std::filesystem::exists(dir)))
        std::filesystem::create_directories(dir);

    filePath << "/logfile_" << calcDateAndTime() << "_" << kernelName << "_vis_" << viscosity << ".txt";
    return filePath.str();
}
