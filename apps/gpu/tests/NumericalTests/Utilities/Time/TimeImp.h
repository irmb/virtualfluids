#ifndef TIME_TRACKING_IMP_H
#define TIME_TRACKING_IMP_H

#include "TimeInfo.h"
#include "TimeTracking.h"

#include <memory>
#include <ctime>

class TimeImp : public TimeTracking, public TimeInfo
{
public:
	static std::shared_ptr<TimeImp> getNewInstance();

	void setSimulationStartTime() override;
	void setSimulationEndTime() override;
	void setTestStartTime() override;
	void setTestEndTime() override;
	void setAnalyticalResultWriteStartTime() override;
	void setAnalyticalResultWriteEndTime() override;
	void setResultCheckStartTime() override;
	void setResultCheckEndTime() override;

	std::string getSimulationTime() override;
	std::string getResultCheckTime() override;
	std::string getTestTime() override;
	std::string getAnalyticalResultWriteTime() override;

private:
	TimeImp();
	double calcSimulationTime();
	float calcResultCheckTime();
	float calcTestTime();
	double calcAnalyticalResultWriteTime();

	time_t simulationStartTime, simulationEndTime;
	clock_t resultCheckStartTime, resultCheckEndTime;
	clock_t testStartTime, testEndTime;
	time_t analyticalResultWriteStartTime, analyticalResultWriteEndTime;
};
#endif