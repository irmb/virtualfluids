#ifndef SIMULATION_LOGFILE_INFORMATION_H
#define SIMULATION_LOGFILE_INFORMATION_H

#include "Utilities\LogFileInformation\LogFileInformationImp.h"

#include <string>

class SimulationLogFileInformation : public LogFileInformationImp
{
public:
	virtual std::string getFilePathExtensionOne() = 0;
	virtual std::string getFilePathExtensionTwo() = 0;

private:

};
#endif 