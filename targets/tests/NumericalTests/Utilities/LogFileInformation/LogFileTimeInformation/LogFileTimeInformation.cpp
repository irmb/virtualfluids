#include "LogFileTimeInformation.h"

#include "Utilities\SimulationInfo\SimulationInfo.h"

#include <iomanip>

std::shared_ptr<LogFileTimeInformation> LogFileTimeInformation::getNewInstance(std::vector<std::shared_ptr<SimulationInfo> > simInfo, bool fileWriting)
{
	return std::shared_ptr<LogFileTimeInformation>(new LogFileTimeInformation(simInfo, fileWriting));
}

std::string LogFileTimeInformation::getOutput()
{
	makeCenterHead("Simulation Time Information");
	oss << "FileWriting=" << std::boolalpha << fileWriting <<std::endl;
	oss << std::endl;
	oss << std::left << std::setfill(' ') << std::setw(11) << "" << "TestName \t \t \t" << " L\t\t" << "Time for Test" << std::endl;
	oss << std::endl;
	for (int i = 0; i < simInfo.size(); i++)
		oss << simInfo.at(i)->getRunTimeOutput();
	oss << std::endl;
	return oss.str();
}

LogFileTimeInformation::LogFileTimeInformation(std::vector<std::shared_ptr<SimulationInfo>> simInfo, bool fileWriting) : simInfo(simInfo), fileWriting(fileWriting)
{

}