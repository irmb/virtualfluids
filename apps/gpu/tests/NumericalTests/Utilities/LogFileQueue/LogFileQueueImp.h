#ifndef LOGFILE_QUEUE_IMP_H
#define LOGFILE_QUEUE_IMP_H

#include "LogFileQueue.h"

#include <string>
#include <vector>
#include <memory>

class LogFileWriter;

class LogFileQueueImp : public LogFileQueue 
{
public:
    static std::shared_ptr<LogFileQueueImp> getNewInstance(std::string basicLogFilePath);

    void writeLogFiles() override;
    void addLogFileWriter(std::shared_ptr<LogFileWriter> aLogFileWriter);

private:
    LogFileQueueImp() = default;
    LogFileQueueImp(std::string basicLogFilePath);

    std::string calcDateAndTime();

    std::string basicLogFilePath;
    std::vector<std::shared_ptr<LogFileWriter> > logFileWriter;
    time_t now;
    struct tm nowLocal;
};
#endif