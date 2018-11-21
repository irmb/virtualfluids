#include "ColorConsoleOutputImp.h"

#include <iomanip>
#include <ctime>

#include "Utilities\SimulationInfo\SimulationInfo.h"


std::shared_ptr<ColorConsoleOutput> ColorConsoleOutputImp::getNewInstance()
{
	return std::shared_ptr<ColorConsoleOutput>(new ColorConsoleOutputImp());
}

void ColorConsoleOutputImp::makeTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3)
{
	setColor(testPassed);
	printTestStart();

	std::ostringstream oss;
	oss << "Kernel: " << simInfo1->getKernelName();
	print(oss.str());
	oss.str(std::string());

	oss << "Viscosity: " << simInfo1->getViscosity();
	print(oss.str());
	oss.str(std::string());

	oss << simInfo1->getSimulationName();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << simInfo1->getLx() << simInfo1->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << simInfo2->getLx() << simInfo2->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	oss << nameWerte1 << ": " << testWert1 << std::setw(5) << "\t" << nameWerte2 << ": " << testWert2;
	print(oss.str());
	oss.str(std::string());

	oss << nameWerte3 << ": " << testWert3;
	print(oss.str());
	oss.str(std::string());

	printTestEnd(testPassed);
}

void ColorConsoleOutputImp::makeSimulationHeadOutput(std::shared_ptr< SimulationInfo> simInfo)
{
	std::ostringstream ossLine1;
	ossLine1 << "# Kernel: " << std::setfill(' ') << std::left << std::setw(38) << simInfo->getKernelName() << "#";

	std::ostringstream ossLine2;
	ossLine2 << "# Viscosity: " << std::setfill(' ') << std::left << std::setw(35) << simInfo->getViscosity() << "#";

	std::ostringstream ossLine3;
	ossLine3 << "# SIMULATION: " << std::setfill(' ') << std::left << std::setw(34) << simInfo->getSimulationName() << "#";
	
	std::ostringstream ossLine4;
	ossLine4 << std::setfill(' ') << std::left << std::setw(14) << "#" << std::setw(34) << simInfo->getSimulationParameterString() << "#";

	std::ostringstream ossLine5;
	ossLine5 << "# L: " << std::setfill(' ') << std::left << std::setw(43) << simInfo->getLx() << "#";

	std::ostringstream ossLine6;
	time_t now;
	struct tm nowLocal;
	now = time(NULL);
	nowLocal = *localtime(&now);
	ossLine6 << "# DATE: " << std::setfill('0') << std::setw(2) << nowLocal.tm_mday << "." << std::setw(2) << nowLocal.tm_mon + 1 << "." << nowLocal.tm_year + 1900 << "   TIME: " << std::setw(2) << nowLocal.tm_hour << ":" << std::setw(2) << nowLocal.tm_min << ":" << std::setw(2) << nowLocal.tm_sec << "\t" << "\t#";

	printGreenHashLine();
	printGreen(ossLine1.str());
	printGreen(ossLine2.str());
	printGreen(ossLine3.str());
	printGreen(ossLine4.str());
	printGreen(ossLine5.str());
	printGreen(ossLine6.str());
	printGreenHashLine();
}

void ColorConsoleOutputImp::makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests)
{
	setColor(numberOfPassedTests == numberOfTests);
	printTestPassed(numberOfPassedTests, numberOfTests);
	printLine();
}

void ColorConsoleOutputImp::makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests)
{
	setColor(numberOfPassedTests == numberOfTests);
	printLine();
	printTestPassed(numberOfPassedTests, numberOfTests);
	std::cout << std::endl;
}

void ColorConsoleOutputImp::printTestStart()
{
	testing::internal::ColoredPrintf(color, "[----------]");
	std::cout << std::endl;
	testing::internal::ColoredPrintf(color, "[Run Test  ]");
	std::cout << std::endl;
	testing::internal::ColoredPrintf(color, "[TestInfo  ]");
	std::cout << std::endl;
}

void ColorConsoleOutputImp::printTestEnd(bool testPassed)
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

void ColorConsoleOutputImp::print(std::string output)
{
	testing::internal::ColoredPrintf(color, "[          ] ");
	testing::internal::ColoredPrintf(testing::internal::COLOR_DEFAULT, output.c_str());
	std::cout << std::endl;
}

void ColorConsoleOutputImp::setColor(bool testPassed)
{
	if (testPassed)
		color = testing::internal::COLOR_GREEN;
	else
		color = testing::internal::COLOR_RED;
}

void ColorConsoleOutputImp::printTestPassed(int numberOfPassedTests, int numberOfTests)
{
	std::ostringstream test;
	test << "[----------]" << std::endl;
	test << "[----------] Test Summary" << std::endl;
	test << "[----------] " << numberOfPassedTests << " out of " << numberOfTests << " tests passed" << std::endl;
	test << "[----------]" << std::endl;
	testing::internal::ColoredPrintf(color, test.str().c_str());
}

void ColorConsoleOutputImp::printLine()
{
	testing::internal::ColoredPrintf(color, "----------------------------------------------------------------------\n");
}

void ColorConsoleOutputImp::printGreen(std::string output)
{
	testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, output.c_str());
	std::cout << std::endl;
}

void ColorConsoleOutputImp::printGreenHashLine()
{
	testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "#################################################\n"); 
}
