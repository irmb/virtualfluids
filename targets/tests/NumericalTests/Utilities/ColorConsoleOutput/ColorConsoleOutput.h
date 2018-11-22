#ifndef COLOR_CONSOLE_OUTPUT_H
#define COLOR_CONSOLE_OUTPUT_H

#include <memory>
#include <sstream>

class SimulationInfo;

class ColorConsoleOutput
{
public:
	virtual void makePhiAndNuTestOutput(bool testPassed, std::shared_ptr< SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3) = 0;
	virtual void makeL2NormTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo, unsigned int basicTimeStep, unsigned int divergentTimeStep, double testWert1, double testWert2, double testWert3) = 0;
	virtual void makeSimulationHeadOutput(std::shared_ptr< SimulationInfo> simInfo) = 0;
	virtual void makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests) = 0;
	virtual void makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests) = 0;

private:

};
#endif