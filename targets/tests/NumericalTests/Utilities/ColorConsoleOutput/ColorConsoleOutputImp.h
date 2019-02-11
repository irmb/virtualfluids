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
		enum GTestColor {
			COLOR_DEFAULT,
			COLOR_RED,
			COLOR_GREEN,
			COLOR_YELLOW
		};

		// in case of unresoved external while using shared libraries
		// add in gtest.h line 167 and 168:
		// enum GTestColor;
		//       void GTEST_API_ ColoredPrintf(GTestColor color, const char* fmt, ...);
		// see commit: 4c0ed885ceab18b9df7a2495c77a51e236aee6f1
		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
	}
}

class ColorConsoleOutputImp : public ColorConsoleOutput
{
public:
	static std::shared_ptr<ColorConsoleOutput> getInstance();

	void makeNyTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string dataToCalc, double nu1, double nu2, double nuDiff1, double nuDiff2, double ooa);
	void makePhiTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string dataToCalc, double phiDiff1, double phiDiff2, double ooa);
	void makeL2NormTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo, unsigned int basicTimeStep, unsigned int divergentTimeStep, std::string dataToCalc, double testWert1, double testWert2, double testWert3);
	void makeL2NormTestErrorOutput(std::string errorMessage, std::shared_ptr<SimulationInfo> simInfo, unsigned int basicTimeStep, unsigned int divergentTimeStep, std::string dataToCalc);
	void makeL2NormBetweenKernelsTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> basicSimInfo, std::shared_ptr<SimulationInfo> divergentSimInfo, std::string dataToCalc, unsigned int timeStep, double l2NormBasicKernel, double l2NormDivergentKernel, double l2NormBetweenKernel);
	void makeSimulationHeadOutput(std::shared_ptr<SimulationInfo> simInfo);
	void makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests);
	void makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests);

private:
	ColorConsoleOutputImp() {};
	void printTestStart();
	void printTestEnd(bool testPassed);
	void printTestEndError();
	void print(std::string output);
	void printColor(std::string output);
	void setColor(bool testPassed);
	void printTestPassed(int numberOfPassedTests, int numberOfTests);
	void printLine();

	void printGreen(std::string output);
	void printGreenHashLine();

	testing::internal::GTestColor color;
};
#endif 