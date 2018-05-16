#include "SimulationInfoImp.h"

#include "Utilities\TestCout\TestCout.h"

#include <sstream>
#include <iomanip>

std::shared_ptr<SimulationInfo> SimulationInfoImp::getNewInstance(std::shared_ptr<TestCout> output, std::string simName, int l)
{
	return std::shared_ptr<SimulationInfo>(new SimulationInfoImp(output, simName, l));
}

void SimulationInfoImp::makeSimulationHeadOutput()
{
	output->makeSimulationHeadOutput(simName, l);
}

void SimulationInfoImp::setStartTime()
{
	startTime = time(NULL);
}

void SimulationInfoImp::setEndTime()
{
	endTime = time(NULL);
}

std::string SimulationInfoImp::getSimulationRunTimeOutput()
{
	std::ostringstream oss;
	oss << std::left << std::setfill(' ') << std::setw(17) << simName << "\t" << std::right << std::setw(3) << l << "\t\t" << std::setw(9) << getSimTime() << " sec" << std::endl;
	return oss.str();
}

SimulationInfoImp::SimulationInfoImp(std::shared_ptr<TestCout> output, std::string simName, int l):simName(simName), l(l), output(output)
{
}

double SimulationInfoImp::getSimTime()
{
	return difftime(endTime, startTime);;
}
