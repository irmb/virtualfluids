#ifndef TEST_COUT_IMP_H
#define TEST_COUT_IMP_H

#include "TestCout.h"

#include <sstream>
#include <memory>
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

class TestCoutImp : public TestCout
{
public:
	static std::shared_ptr<TestCout> getNewInstance();

	void makeTestOutput(bool testPassed, std::string testName, int l1, int l2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3);
	void makeSimulationHeadOutput(std::string simName, int l);
	void makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests);
	void makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests);

private:
	TestCoutImp() {};
	void printTestStart();
	void printTestEnd(bool testPassed);
	void print(std::string output);
	void setColor(bool testPassed);
	void printTestPassed(int numberOfPassedTests, int numberOfTests);
	void printLine();

	void printGreen(std::string output);
	void printGreenHashLine();

	testing::internal::GTestColor color;
};
#endif 