#ifndef SIMULATION_LOGFILE_INFORMATION_H
#define SIMULATION_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"

#include <string>
#include <vector>

class SimulationLogFileInformation : public LogFileInformationImp
{
public:
    virtual ~SimulationLogFileInformation() = default;
    virtual std::string getOutput() = 0;

    virtual std::vector<std::string> getFilePathExtension() = 0;

private:

};
#endif 