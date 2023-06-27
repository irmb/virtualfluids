#ifndef LOG_FILE_WRITER_H
#define LOG_FILE_WRITER_H

#include <string>

class LogFileWriter
{
public:
	virtual ~LogFileWriter() = default;
	virtual void writeLogFile(std::string basicFilePath) = 0;
	
private:

};
#endif 
