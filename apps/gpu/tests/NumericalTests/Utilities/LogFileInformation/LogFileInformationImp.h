#ifndef LOGFILE_INFORMATION_IMP_H
#define LOGFILE_INFORMATION_IMP_H

#include "LogFileInformation.h"

#include <sstream>

class LogFileInformationImp : public LogFileInformation
{
public:
    LogFileInformationImp() = default;
    virtual std::string getOutput() = 0;

protected:
    void makeCenterHead(std::string head);

    std::ostringstream oss;
private:
    void makeHashLine();

};
#endif