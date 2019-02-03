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
	oss << "VTKFileWriting=" << std::boolalpha << fileWriting <<std::endl << std::endl;
	for (int i = 0; i < simInfo.size(); i++) {
		oss << simInfo.at(i)->getRunTimeOutput() << std::endl;
	}
	return oss.str();
}

LogFileTimeInformation::LogFileTimeInformation(std::vector<std::shared_ptr<SimulationInfo>> simInfo, bool fileWriting) : simInfo(simInfo), fileWriting(fileWriting)
{

}