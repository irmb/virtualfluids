#include "SimulationInfoImp.h"

#include "Utilities/Time/TimeInfo.h"

#include <iomanip>
#include <sstream>

void SimulationInfoImp::setTimeInfo(std::shared_ptr<TimeInfo> timeInfo)
{
	this->timeInfo = timeInfo;
}

std::string SimulationInfoImp::getKernelName()
{
	return kernelName;
}

double SimulationInfoImp::getViscosity()
{
	return viscosity;
}

std::string SimulationInfoImp::getSimulationName()
{
	return simulationName;
}

std::string SimulationInfoImp::getSimulationParameterString()
{
	return simulationParameterString;
}

int SimulationInfoImp::getLx()
{
	return lx;
}

int SimulationInfoImp::getNumberOfSimulations()
{
	return numberOfSimulations;
}

int SimulationInfoImp::getSimulationID()
{
	return simID;
}

std::string SimulationInfoImp::getRunTimeOutput()
{
	std::ostringstream oss;
	oss << "SimulationTime_" << lx << "=" << timeInfo->getSimulationTime() << std::endl;
	oss << "TestTime_" << lx << "=" << timeInfo->getTestTime() << std::endl;
	oss << "AnalyticalVTKFileWritingTime_" << lx << "=" << timeInfo->getAnalyticalResultWriteTime() << std::endl;
	return oss.str();
}

SimulationInfoImp::SimulationInfoImp(int simID, std::string kernelName, double viscosity, int lx, int numberOfSimulations, std::string simulationName)
	: simID(simID), lx(lx), viscosity(viscosity), kernelName(kernelName), numberOfSimulations(numberOfSimulations), simulationName(simulationName)
{
}
