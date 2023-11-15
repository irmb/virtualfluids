#ifndef COLOR_CONSOLE_OUTPUT_IMP_H
#define COLOR_CONSOLE_OUTPUT_IMP_H

#include "ColorConsoleOutput.h"

#include <iostream>
#include <string>


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

    // not used at the moment
    std::string color;
};
#endif