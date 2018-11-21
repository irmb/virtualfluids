#ifndef TEST_COUT_H
#define TEST_COUT_H

#include <iostream>
#include <memory>

class SimulationInfo;

class TestCout
{
public:
	virtual void makeTestOutput(bool testPassed, std::shared_ptr< SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3) = 0;
	virtual void makeSimulationHeadOutput(std::shared_ptr< SimulationInfo> simInfo) = 0;
	virtual void makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests) = 0;
	virtual void makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests) = 0;
private:
};
#endif 