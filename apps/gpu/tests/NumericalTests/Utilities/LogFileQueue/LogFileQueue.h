#ifndef LOGFILE_QUEUE_H
#define LOGFILE_QUEUE_H

#include <memory>

class LogFileWriter;

class LogFileQueue
{
public:
	virtual void writeLogFiles() = 0;

private:

};
#endif