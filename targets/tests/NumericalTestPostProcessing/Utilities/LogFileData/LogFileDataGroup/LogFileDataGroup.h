#ifndef LOG_FILE_DATA_GROUP_H
#define LOG_FILE_DATA_GROUP_H

#include <memory>

class LogFileData;

class LogFileDataGroup
{
public:
	virtual std::shared_ptr<LogFileData> getLogFileData(int number) = 0;
	virtual int getGroupSize() = 0;
};
#endif