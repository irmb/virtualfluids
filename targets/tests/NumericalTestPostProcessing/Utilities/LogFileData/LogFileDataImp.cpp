#include "LogFileDataImp.h"

std::shared_ptr<LogFileDataImp> LogFileDataImp::getNewInstance()
{
	return std::shared_ptr<LogFileDataImp>(new LogFileDataImp());
}

void LogFileDataImp::setBasisTimeStepLength(int basisTimeStepLength)
{
	this->basisTimeStepLength = basisTimeStepLength;
}

void LogFileDataImp::setKernel(std::string kernelName)
{
	this->kernelName = kernelName;
}

void LogFileDataImp::setNumberOfTimeSteps(int numberOfTimeSteps)
{
	this->numberOfTimeSteps = numberOfTimeSteps;
}

void LogFileDataImp::setViscosity(double viscosity)
{
	this->viscosity = viscosity;
}

int LogFileDataImp::getBasisTimeStepLength()
{
	return basisTimeStepLength;
}

std::string LogFileDataImp::getKernel()
{
	return kernelName;
}

int LogFileDataImp::getNumberOfTimeSteps()
{
	return numberOfTimeSteps;
}

double LogFileDataImp::getViscosity()
{
	return viscosity;
}

LogFileDataImp::LogFileDataImp()
{
}
