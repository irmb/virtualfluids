#ifndef TIME_TRACKING_IMP_H
#define TIME_TRACKING_IMP_H

#include "TimeInfo.h"
#include "TimeTracking.h"

#include <memory>
#include <time.h>

class TimeImp : public TimeTracking, public TimeInfo
{
public:
	static std::shared_ptr<TimeImp> getNewInstance();

	void setSimulationStartTime();
	void setSimulationEndTime();
	void setTestStartTime();
	void setTestEndTime();
	void setAnalyticalResultWriteStartTime();
	void setAnalyticalResultWriteEndTime();

	std::string getSimulationTime();
	std::string getTestTime();
	std::string getAnalyticalResultWriteTime();

private:
	TimeImp();
	double calcSimulationTime();
	float calcTestTime();
	double calcAnalyticalResultWriteTime();

	time_t simulationStartTime, simulationEndTime;
	clock_t testStartTime, testEndTime;
	time_t analyticalResultWriteStartTime, analyticalResultWriteEndTime;
};
#endif