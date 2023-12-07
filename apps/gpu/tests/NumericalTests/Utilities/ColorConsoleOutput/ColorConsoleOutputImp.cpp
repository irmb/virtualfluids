#include "ColorConsoleOutputImp.h"

#include "Utilities/SimulationInfo/SimulationInfo.h"

#include <ctime>
#include <iomanip>

#include <logger/Logger.h>

#include <basics/DataTypes.h>

void log(const char* fmt) {
    VF_LOG_INFO("{}", fmt);
}

std::shared_ptr<ColorConsoleOutput> ColorConsoleOutputImp::getInstance()
{
    static std::shared_ptr<ColorConsoleOutput> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<ColorConsoleOutput>(new ColorConsoleOutputImp());
    return uniqueInstance;
}

void ColorConsoleOutputImp::makeSimulationHeadOutput(std::shared_ptr<SimulationInfo> simInfo)
{
    std::ostringstream ossLine0;
    ossLine0 << "# Simulation Number " << simInfo->getSimulationID() << " of " << simInfo->getNumberOfSimulations();
    int length = 49 - ossLine0.str().size();
    ossLine0 << std::setfill(' ') << std::right << std::setw(length) << "#";

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

void ColorConsoleOutputImp::makeTestOutput(std::vector<std::string> testOutput, TestStatus status)
{
    setColor(status);
    printTestStart();

    printColor("");
    printColor(testOutput.at(0));
    printColor("");

    for (uint i = 1; i < testOutput.size(); i++)
        print(testOutput.at(i));

    printColor("");
    printColor(testOutput.at(0));
    printColor("");

    printTestEnd(status);
}

void ColorConsoleOutputImp::makeFinalTestOutputHead(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest)
{
    setColor(numberOfTests == numberOfPassedTest);
    printTestPassed(numberOfTests, numberOfExecutedTest, numberOfPassedTest, numberOfFailedTest, numberOfErrorTest, numberOfNotExecutedTest);
    printLine();
}

void ColorConsoleOutputImp::makeFinalTestOutputFoot(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest)
{
    setColor(numberOfTests == numberOfPassedTest);
    printLine();
    printTestPassed(numberOfTests, numberOfExecutedTest, numberOfPassedTest, numberOfFailedTest, numberOfErrorTest, numberOfNotExecutedTest);
}

void ColorConsoleOutputImp::printTestStart()
{
    log( "[-----------]");
    log( "[Run Test   ]");
    log( "[TestInfo   ]");
}

void ColorConsoleOutputImp::printTestEnd(TestStatus status)
{
    log( "[   TestInfo]");
    switch (status)
    {
    case passed: log( "[     PASSED]");
        break;
    case failed: log( "[     FAILED]");
        break;
    case test_error: log( "[      ERROR]");
        break;
    case simulationCrashed: log( "[Sim crashed]");
        break;
    default:
        break;
    }

    log( "[-----------]");
}

void ColorConsoleOutputImp::print(std::string output)
{
    log("[           ] ");
    log(output.c_str());
}

void ColorConsoleOutputImp::printColor(std::string output)
{
    log("[-----------] ");
    log(output.c_str());
}

void ColorConsoleOutputImp::setColor(TestStatus status)
{
    switch (status)
    {
    case passed: color = "green";
        break;
    case failed: color = "red";
        break;
    case test_error: color = "yellow";
        break;
    case simulationCrashed: color = "yellow";
        break;
    default:
        break;
    }
}

void ColorConsoleOutputImp::setColor(bool passed)
{
    if (passed)
        color = "green";
    else
        color = "red";
}

void ColorConsoleOutputImp::printTestPassed(int numberOfTests, int numberOfExecutedTest, int numberOfPassedTest, int numberOfFailedTest, int numberOfErrorTest, int numberOfNotExecutedTest)
{
    std::ostringstream test;
    test << "[-----------]" << std::endl;
    test << "[-----------] Test Summary" << std::endl;
    test << "[-----------] " << numberOfTests << " initialized Tests" << std::endl;
    test << "[-----------]" << std::endl;
    test << "[-----------] " << numberOfExecutedTest << " out of " << numberOfTests << " Tests executed" << std::endl;
    test << "[-----------] " << numberOfErrorTest << " out of " << numberOfTests << " Tests executed and completed with error" << std::endl;
    test << "[-----------] " << numberOfNotExecutedTest << " out of " << numberOfTests << " Tests not executed" << std::endl;
    test << "[-----------]" << std::endl;
    test << "[-----------] " << numberOfPassedTest << " out of " << numberOfExecutedTest << " executed Tests passed" << std::endl;
    test << "[-----------] " << numberOfFailedTest << " out of " << numberOfExecutedTest << " executed Tests failed" << std::endl;
    test << "[-----------]" << std::endl;
    log(test.str().c_str());
}

void ColorConsoleOutputImp::printLine()
{
    log("----------------------------------------------------------------------");
}

void ColorConsoleOutputImp::printGreen(std::string output)
{
    log(output.c_str());
}

void ColorConsoleOutputImp::printGreenHashLine()
{
    log("#################################################");
}
