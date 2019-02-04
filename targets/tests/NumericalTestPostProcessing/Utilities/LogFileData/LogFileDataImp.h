#ifndef LOGFILE_DATA_IMP_H
#define LOGFILE_DATA_IMP_H

#include "LogFileData.h"

#include <memory>

class LogFileDataImp : public LogFileData
{
public:
	static std::shared_ptr<LogFileDataImp> getNewInstance();

	int getBasisTimeStepLength();
	std::string getKernel();
	int getNumberOfTimeSteps();
	double getViscosity();
	
	void setBasisTimeStepLength(int basisTimeStepLength);
	void setKernel(std::string kernelName);
	void setNumberOfTimeSteps(int numberOfTimeSteps);
	void setViscosity(double viscosity);

private:
	LogFileDataImp();

	int basisTimeStepLength;
	std::string kernelName;
	int numberOfTimeSteps;
	double viscosity;
};
#endif