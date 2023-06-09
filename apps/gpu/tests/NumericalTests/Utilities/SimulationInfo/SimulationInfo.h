#ifndef SIMULATION_INFO_H
#define SIMULATION_INFO_H

#include <memory>
#include <string>
#include <vector>

class TimeInfo;

class SimulationInfo
{
public:
	virtual ~SimulationInfo() = default;
	virtual std::string getKernelName() = 0;
	virtual double getViscosity() = 0;
	virtual std::string getSimulationName() = 0;
	virtual std::string getSimulationParameterString() = 0;
	virtual int getLx() = 0;
	virtual int getNumberOfSimulations() = 0;
	virtual int getSimulationID() = 0;
	virtual std::string getRunTimeOutput() = 0;
	virtual std::vector<std::string> getDataToCalcTests() = 0; 

	virtual void setTimeInfo(std::shared_ptr<TimeInfo> timeInfo) = 0;
};
#endif 