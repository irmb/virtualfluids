#ifndef LOG_FILE_DATA_GROUP_IMP_H
#define LOG_FILE_DATA_GROUP_IMP_H

#include "LogFileDataGroup.h"

#include <vector>

class LogFileDataGroupImp : public LogFileDataGroup
{
public:
    static std::shared_ptr<LogFileDataGroupImp> getNewInstance();

    std::shared_ptr<LogFileData> getLogFileData(int number);
    int getGroupSize();
    void addLogFileData(std::shared_ptr<LogFileData> logFileData);

private:
    LogFileDataGroupImp();

    std::vector<std::shared_ptr<LogFileData>> data;
};
#endif