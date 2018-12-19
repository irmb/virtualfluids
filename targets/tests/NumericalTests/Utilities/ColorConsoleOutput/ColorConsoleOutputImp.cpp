#include "ColorConsoleOutputImp.h"

#include <iomanip>
#include <ctime>

#include "Utilities\SimulationInfo\SimulationInfo.h"


std::shared_ptr<ColorConsoleOutput> ColorConsoleOutputImp::getInstance()
{
	static std::shared_ptr<ColorConsoleOutput> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<ColorConsoleOutput>(new ColorConsoleOutputImp());
	return uniqueInstance;
}

void ColorConsoleOutputImp::makeNuTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string dataToCalc, double nu1, double nu2, double nuDiff1, double nuDiff2, double ooa)
{
	setColor(testPassed);
	printTestStart();

	printColor("");
	printColor("Nu Test");
	printColor("");

	std::ostringstream oss;
	oss << "Kernel: " << simInfo1->getKernelName();
	print(oss.str());
	oss.str(std::string());

	oss << "Viscosity: " << simInfo1->getViscosity();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << simInfo1->getSimulationName();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << std::setfill(' ') << std::right << std::setw(4) << simInfo1->getLx() << simInfo1->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << std::setfill(' ') << std::right << std::setw(4) << simInfo2->getLx() << simInfo2->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "DataToCalculate: " << dataToCalc;
	print(oss.str());
	oss.str(std::string());

	oss << "StartTimeStep: " << startTimeStep;
	print(oss.str());
	oss.str(std::string());

	oss << "EndTimeStep: " << endTimeStep;
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "Nu" << simInfo1->getLx() << ": " << nu1;
	print(oss.str());
	oss.str(std::string());
	oss << "Nu" << simInfo2->getLx() << ": " << nu2;
	print(oss.str());
	oss.str(std::string());
	oss << "NuDiff" << simInfo1->getLx() << ": " << nuDiff1;
	print(oss.str());
	oss.str(std::string());
	oss << "NuDiff" << simInfo2->getLx() << ": " << nuDiff2;
	print(oss.str());
	oss.str(std::string());
	oss << "OrderOfAccuracy: " << ooa;
	print(oss.str());
	oss.str(std::string());

	printColor("");
	printColor("Nu Test");
	printColor("");

	printTestEnd(testPassed);
}

void ColorConsoleOutputImp::makePhiTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string dataToCalc, double phiDiff1, double phiDiff2, double ooa)
{
	setColor(testPassed);
	printTestStart();

	printColor("");
	printColor("Phi Test");
	printColor("");

	std::ostringstream oss;
	oss << "Kernel: " << simInfo1->getKernelName();
	print(oss.str());
	oss.str(std::string());

	oss << "Viscosity: " << simInfo1->getViscosity();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << simInfo1->getSimulationName();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << std::setfill(' ') << std::right << std::setw(4) << simInfo1->getLx() << simInfo1->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << std::setfill(' ') << std::right << std::setw(4) << simInfo2->getLx() << simInfo2->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "DataToCalculate: " << dataToCalc;
	print(oss.str());
	oss.str(std::string());

	oss << "StartTimeStep: " << startTimeStep;
	print(oss.str());
	oss.str(std::string());

	oss << "EndTimeStep: " << endTimeStep;
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "PhiDiff" << simInfo1->getLx() << ": " << phiDiff1;
	print(oss.str());
	oss.str(std::string());
	oss << "PhiDiff" << simInfo2->getLx() << ": " << phiDiff2;
	print(oss.str());
	oss.str(std::string());
	oss << "OrderOfAccuracy: " << ooa;
	print(oss.str());
	oss.str(std::string());

	printColor("");
	printColor("Phi Test");
	printColor("");

	printTestEnd(testPassed);
}

void ColorConsoleOutputImp::makeL2NormTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo, unsigned int basicTimeStep, unsigned int divergentTimeStep, std::string dataToCalc, double testWert1, double testWert2, double testWert3)
{
	setColor(testPassed);
	printTestStart();

	printColor("");
	printColor("L2 Norm Test");
	printColor("");

	std::ostringstream oss;
	oss << "Kernel: " << simInfo->getKernelName();
	print(oss.str());
	oss.str(std::string());

	oss << "Viscosity: " << simInfo->getViscosity();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << simInfo->getSimulationName();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << simInfo->getLx() << simInfo->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "DataToCalculate: " << dataToCalc;
	print(oss.str());
	oss.str(std::string());
	oss << "BasicTimeStep: " << basicTimeStep;
	print(oss.str());
	oss.str(std::string());

	oss << "DivergentTimeStep: " << divergentTimeStep;
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "L2Norm BasicTimeStep: " << testWert1;
	print(oss.str());
	oss.str(std::string());

	oss << "L2Norm DivergentTimeStep: " << testWert2;
	print(oss.str());
	oss.str(std::string());

	oss << "L2NormDiff: " << testWert3;
	print(oss.str());
	oss.str(std::string());

	printColor("");
	printColor("L2 Norm Test");
	printColor("");

	printTestEnd(testPassed);
}

void ColorConsoleOutputImp::makeL2NormBetweenKernelsTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> basicSimInfo, std::shared_ptr<SimulationInfo> divergentSimInfo, std::string dataToCalc, unsigned int timeStep, double l2NormBasicKernel, double l2NormDivergentKernel, double l2NormBetweenKernel)
{
	setColor(testPassed);
	printTestStart();

	printColor("");
	printColor("L2 Norm Between Kernels Test");
	printColor("");

	std::ostringstream oss;
	oss << "Basic Kernel: " << basicSimInfo->getKernelName();
	print(oss.str());
	oss.str(std::string());
	oss << "Divergent Kernel: " << divergentSimInfo->getKernelName();
	print(oss.str());
	oss.str(std::string());

	oss << "Viscosity: " << basicSimInfo->getViscosity();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << basicSimInfo->getSimulationName();
	print(oss.str());
	oss.str(std::string());

	oss << "L: " << basicSimInfo->getLx() << basicSimInfo->getSimulationParameterString();
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "DataToCalculate: " << dataToCalc;
	print(oss.str());
	oss.str(std::string());
	oss << "TimeStep: " << timeStep;
	print(oss.str());
	oss.str(std::string());

	print(oss.str());
	oss << "L2Norm BasicKernel: " << l2NormBasicKernel;
	print(oss.str());
	oss.str(std::string());

	oss << "L2Norm DivergentKernel: " << l2NormDivergentKernel;
	print(oss.str());
	oss.str(std::string());

	oss << "L2NormDiff: " << l2NormBetweenKernel;
	print(oss.str());
	oss.str(std::string());

	printColor("");
	printColor("L2 Norm Between Kernels Test");
	printColor("");

	printTestEnd(testPassed);
}

void ColorConsoleOutputImp::makeSimulationHeadOutput(std::shared_ptr< SimulationInfo> simInfo)
{
	std::ostringstream ossLine0;
	ossLine0 << "# Simulation Number " << simInfo->getSimulationID() << " of " << std::setfill(' ') << std::left << std::setw(23) << simInfo->getNumberOfSimulations() << "#";

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
	printGreen(ossLine0.str());
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

void ColorConsoleOutputImp::printColor(std::string output)
{
	testing::internal::ColoredPrintf(color, "[----------] ");
	testing::internal::ColoredPrintf(color, output.c_str());
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
