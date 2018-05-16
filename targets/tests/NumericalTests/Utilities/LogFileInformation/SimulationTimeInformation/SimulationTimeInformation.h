#ifndef SIMULATION_TIME_INFORMATION_H
#define SIMULATION_TIME_INFORMATION_H

#include "../LogFileInformationImp.h"

#include <memory>
#include <vector>

class SimulationInfo;

class SimulationTimeInformation : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileInformation> getNewInstance(std::vector< std::shared_ptr< SimulationInfo> > simInfo, bool fileWriting);
	std::string getOutput();

private:
	SimulationTimeInformation();
	SimulationTimeInformation(std::vector< std::shared_ptr< SimulationInfo> > simInfo, bool fileWriting);

	std::vector< std::shared_ptr< SimulationInfo> > simInfo;
	bool fileWriting;
};
#endif