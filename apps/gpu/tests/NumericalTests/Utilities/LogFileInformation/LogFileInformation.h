#ifndef LOGFILE_INFORMATION_H
#define LOGFILE_INFORMATION_H

#include <string>

class LogFileInformation
{
public:
    virtual ~LogFileInformation() = default;
    virtual std::string getOutput() = 0;

private:

};
#endif