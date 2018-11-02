#ifndef TEST_COUT_IMP_H
#define TEST_COUT_IMP_H

#include "TestCout.h"
#include "ColorOutput.h"

#include <sstream>
#include <memory>


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