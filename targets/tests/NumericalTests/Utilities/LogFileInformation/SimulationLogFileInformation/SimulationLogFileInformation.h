#ifndef SIMULATION_LOGFILE_INFORMATION_H
#define SIMULATION_LOGFILE_INFORMATION_H

#include "Utilities\LogFileInformation\LogFileInformation.h"

#include <string>

class SimulationLogFileInformation : public LogFileInformation
{
public:
	virtual std::string getFilePathExtension() = 0;

private:

};
#endif 