#ifndef COLOR_CONSOLE_OUTPUT_IMP_H
#define COLOR_CONSOLE_OUTPUT_IMP_H

#include "ColorConsoleOutput.h"

#include <iostream>

#include <gtest/gtest.h>

// https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework/29155677#29155677
namespace testing
{
	namespace internal
	{
		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
	}
}

class ColorConsoleOutputImp : public ColorConsoleOutput
{
public:
	static std::shared_ptr<ColorConsoleOutput> getInstance();

	void makeSimulationHeadOutput(std::shared_ptr<SimulationInfo> simInfo);
	void makeTestOutput(std::vector<std::string> testOutput, TestStatus status);
	void makeFinalTestOutputHead(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest);
	void makeFinalTestOutputFoot(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest);

private:
	ColorConsoleOutputImp() {};
	void printTestStart();
	void printTestEnd(TestStatus status);
	void print(std::string output);
	void printColor(std::string output);
	void setColor(TestStatus status);
	void setColor(bool passed);
	void printTestPassed(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest);
	void printLine();

	void printGreen(std::string output);
	void printGreenHashLine();

	testing::internal::GTestColor color;
};
#endif 