#include "TestCoutImp.h"

#include <iomanip>
#include <ctime>


std::shared_ptr<TestCout> TestCoutImp::getNewInstance()
{
	return std::shared_ptr<TestCout>(new TestCoutImp());
}

void TestCoutImp::makeTestOutput(bool testPassed, std::string testName, int l1, int l2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3)
{
	setColor(testPassed);
	printTestStart();

	std::ostringstream oss1;
	oss1 << testName;
	print(oss1.str());

	std::ostringstream oss2;
	oss2 << "L: " << l1 << "\t" << "\t" << "\t" << "L: " << l2;
	print(oss2.str());

	std::ostringstream oss3;
	oss3 << nameWerte1 << ": " << testWert1 << std::setw(5) << "\t" << nameWerte2 << ": " << testWert2;
	print(oss3.str());

	std::ostringstream oss4;
	oss4 << nameWerte3 << ": " << testWert3;
	print(oss4.str());
	
	printTestEnd(testPassed);
}

void TestCoutImp::makeSimulationHeadOutput(std::string simName, int l)
{
	std::ostringstream oss1;
	oss1 << "# SIMULATION: " << std::setfill(' ') << std::left << std::setw(34) << simName << "#";
	
	std::ostringstream oss2;
	oss2 << "# L: " << std::setfill(' ') << std::left << std::setw(43) << l << "#";

	std::ostringstream oss3;
	time_t now;
	struct tm nowLocal;
	now = time(NULL);
	nowLocal = *localtime(&now);
	oss3 << "# DATE: " << std::setfill('0') << std::setw(2) << nowLocal.tm_mday << "." << std::setw(2) << nowLocal.tm_mon + 1 << "." << nowLocal.tm_year + 1900 << "   TIME: " << std::setw(2) << nowLocal.tm_hour << ":" << std::setw(2) << nowLocal.tm_min << ":" << std::setw(2) << nowLocal.tm_sec << "\t" << "\t#";

	printGreenHashLine();
	printGreen(oss1.str());
	printGreen(oss2.str());
	printGreen(oss3.str());
	printGreenHashLine();
}

void TestCoutImp::makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests)
{
	setColor(numberOfPassedTests == numberOfTests);
	printTestPassed(numberOfPassedTests, numberOfTests);
	printLine();
}

void TestCoutImp::makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests)
{
	setColor(numberOfPassedTests == numberOfTests);
	printLine();
	printTestPassed(numberOfPassedTests, numberOfTests);
	std::cout << std::endl;
}

void TestCoutImp::printTestStart()
{
	testing::internal::ColoredPrintf(color, "[----------]");
	std::cout << std::endl;
	testing::internal::ColoredPrintf(color, "[Run Test  ]");
	std::cout << std::endl;
	testing::internal::ColoredPrintf(color, "[TestInfo  ]");
	std::cout << std::endl;
}

void TestCoutImp::printTestEnd(bool testPassed)
{
	testing::internal::ColoredPrintf(color, "[  TestInfo]");
	std::cout << std::endl;
	if (testPassed)
		testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[    PASSED]");
	else
		testing::internal::ColoredPrintf(testing::internal::COLOR_RED, "[    FAILED]");
	std::cout << std::endl;
	testing::internal::ColoredPrintf(color, "[----------]");
	std::cout << std::endl << std::endl;
}

void TestCoutImp::print(std::string output)
{
	testing::internal::ColoredPrintf(color, "[          ] ");
	testing::internal::ColoredPrintf(testing::internal::COLOR_DEFAULT, output.c_str());
	std::cout << std::endl;
}

void TestCoutImp::setColor(bool testPassed)
{
	if (testPassed)
		color = testing::internal::COLOR_GREEN;
	else
		color = testing::internal::COLOR_RED;
}

void TestCoutImp::printTestPassed(int numberOfPassedTests, int numberOfTests)
{
	std::ostringstream test;
	test << "[----------]" << std::endl;
	test << "[----------] Test Summary" << std::endl;
	test << "[----------] " << numberOfPassedTests << " out of " << numberOfTests << " tests passed" << std::endl;
	test << "[----------]" << std::endl;
	testing::internal::ColoredPrintf(color, test.str().c_str());
}

void TestCoutImp::printLine()
{
	testing::internal::ColoredPrintf(color, "----------------------------------------------------------------------\n");
}

void TestCoutImp::printGreen(std::string output)
{
	testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, output.c_str());
	std::cout << std::endl;
}

void TestCoutImp::printGreenHashLine()
{
	testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "#################################################\n"); 
}
