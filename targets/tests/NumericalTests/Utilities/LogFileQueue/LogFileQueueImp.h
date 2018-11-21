#ifndef LOGFILE_QUEUE_IMP_H
#define LOGFILE_QUEUE_IMP_H

#include "LogFileQueue.h"

#include <string>
#include <vector>

class LogFileQueueImp : public LogFileQueue 
{
public:
	static std::shared_ptr<LogFileQueueImp> getNewInstance(std::string basicLogFilePath);

	void writeLogFiles();
	void addLogFileWriter(std::shared_ptr< LogFileWriter> aLogFileWriter);

private:
	LogFileQueueImp() {};
	LogFileQueueImp(std::string basicLogFilePath);

	std::string calcDateAndTime();

	std::string basicLogFilePath;
	std::vector< std::shared_ptr< LogFileWriter>> logFileWriter;
	time_t now;
	struct tm nowLocal;
};
#endif