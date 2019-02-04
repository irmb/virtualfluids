#ifndef LOGFILE_DATA_H
#define LOGFILE_DATA_H

#include <string>

class LogFileData
{
public:
	virtual int getBasisTimeStepLength() = 0;
	virtual std::string getKernel() = 0;
	virtual int getNumberOfTimeSteps() = 0;
	virtual double getViscosity() = 0;
};
#endif