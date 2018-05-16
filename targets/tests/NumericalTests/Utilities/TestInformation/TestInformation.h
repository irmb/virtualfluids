#ifndef TEST_INFORMATION_H
#define TEST_INFORMATION_H

class TestInformation
{
public:
	virtual void makeSimulationHeadOutput(int i) = 0;
	virtual void setSimulationStartTime(int i) = 0;
	virtual void setSimulationEndTime(int i) = 0;
	virtual void writeLogFile() = 0;
	virtual void makeFinalTestOutput() = 0;

private:

};
#endif 
