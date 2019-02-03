#include "SimulationInfoImp.h"

#include "Utilities\Time\TimeInfo.h"

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
	oss << std::left << std::setfill(' ') << std::setw(11) << "Simulation" << std::setw(17) << simulationName << "\t" << std::right << std::setw(3) << lx << "\t\t" << std::setw(9) << timeInfo->getSimulationTime() << std::endl;
	oss << std::left << std::setfill(' ') << std::setw(11) << "Test" << std::setw(17) << simulationName << "\t" << std::right << std::setw(3) << lx << "\t\t" << std::setw(9) << timeInfo->getTestTime() << std::endl;
	oss << std::left << std::setfill(' ') << std::setw(11) << "FileWriting" << std::setw(17) << simulationName << "\t" << std::right << std::setw(3) << lx << "\t\t" << std::setw(9) << timeInfo->getAnalyticalResultWriteTime() << std::endl;
	return oss.str();
}

SimulationInfoImp::SimulationInfoImp(int simID, std::string kernelName, double viscosity, int lx, int numberOfSimulations, std::string simulationName)
	: simID(simID), lx(lx), viscosity(viscosity), kernelName(kernelName), numberOfSimulations(numberOfSimulations), simulationName(simulationName)
{
}
