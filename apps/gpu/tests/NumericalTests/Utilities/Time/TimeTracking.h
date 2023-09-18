#ifndef TIME_TRACKING_H
#define TIME_TRACKING_H

class TimeTracking
{
public:
	virtual ~TimeTracking() = default;
	virtual void setSimulationStartTime() = 0;
	virtual void setSimulationEndTime() = 0;

	virtual void setResultCheckStartTime() = 0;
	virtual void setResultCheckEndTime() = 0;

	virtual void setTestStartTime() = 0;
	virtual void setTestEndTime() = 0;

	virtual void setAnalyticalResultWriteStartTime() = 0;
	virtual void setAnalyticalResultWriteEndTime() = 0;
};
#endif
