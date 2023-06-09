#ifndef LOGFILE_QUEUE_H
#define LOGFILE_QUEUE_H

class LogFileQueue
{
public:
	virtual ~LogFileQueue() = default;
	virtual void writeLogFiles() = 0;

private:

};
#endif