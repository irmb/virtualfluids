#ifndef TEST_LOGFILE_INFORMATION_H
#define TEST_LOGFILE_INFORMATION_H

#include "../LogFileInformationImp.h" 

#include <iostream>

class TestLogFileInformation : public LogFileInformationImp
{
public:
    virtual std::string getOutput() = 0;

private:

};
#endif