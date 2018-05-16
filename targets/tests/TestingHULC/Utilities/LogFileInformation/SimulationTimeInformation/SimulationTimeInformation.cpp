#include "SimulationTimeInformation.h"

#include "Utilities/SimulationInfo/SimulationInfo.h"

std::shared_ptr<LogFileInformation> SimulationTimeInformation::getNewInstance(std::vector<std::shared_ptr<SimulationInfo>> simInfo, bool fileWriting)
{
	return std::shared_ptr<LogFileInformation>(new SimulationTimeInformation(simInfo, fileWriting));
}

std::string SimulationTimeInformation::getOutput()
{
	makeCenterHead("Simulation Time Information");
	oss << "FileWriting: " << fileWriting <<std::endl;
	oss << std::endl;
	oss << "TestName \t \t \t" << " L\t\t" << "Time for Test" << std::endl;
	oss << std::endl;
	for (int i = 0; i < simInfo.size(); i++)
		oss << simInfo.at(i)->getSimulationRunTimeOutput();
	oss << std::endl;
	return oss.str();
}

SimulationTimeInformation::SimulationTimeInformation(std::vector<std::shared_ptr<SimulationInfo>> simInfo, bool fileWriting) : simInfo(simInfo), fileWriting(fileWriting)
{
}