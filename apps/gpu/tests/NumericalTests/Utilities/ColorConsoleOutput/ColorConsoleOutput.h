#ifndef COLOR_CONSOLE_OUTPUT_H
#define COLOR_CONSOLE_OUTPUT_H

#include "Utilities/Test/TestStatus.h"

#include <memory>
#include <sstream>
#include <vector>

class SimulationInfo;

class ColorConsoleOutput
{
public:
	virtual void makeSimulationHeadOutput(std::shared_ptr<SimulationInfo> simInfo) = 0;
	virtual void makeTestOutput(std::vector<std::string> testOutput, TestStatus status) = 0;
	virtual void makeFinalTestOutputHead(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest) = 0;
	virtual void makeFinalTestOutputFoot(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest) = 0;
};
#endif