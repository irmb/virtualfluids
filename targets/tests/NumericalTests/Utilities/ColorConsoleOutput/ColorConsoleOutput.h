#ifndef COLOR_CONSOLE_OUTPUT_H
#define COLOR_CONSOLE_OUTPUT_H

#include <memory>
#include <sstream>

class SimulationInfo;

class ColorConsoleOutput
{
public:
	virtual void makeNyTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string dataToCalc, double nu1, double nu2, double nuDiff1, double nuDiff2, double ooa) = 0;
	virtual void makePhiTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo1, std::shared_ptr<SimulationInfo> simInfo2, unsigned int startTimeStep, unsigned int endTimeStep, std::string dataToCalc, double phiDiff1, double phiDiff2, double ooa) = 0;
	virtual void makeL2NormTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> simInfo, unsigned int basicTimeStep, unsigned int divergentTimeStep, std::string dataToCalc, double testWert1, double testWert2, double testWert3) = 0;
	virtual void makeL2NormBetweenKernelsTestOutput(bool testPassed, std::shared_ptr<SimulationInfo> basicSimInfo, std::shared_ptr<SimulationInfo> divergentSimInfo, std::string dataToCalc, unsigned int timeStep, double l2NormBasicKernel, double l2NormDivergentKernel, double l2NormBetweenKernel) = 0;
	virtual void makeSimulationHeadOutput(std::shared_ptr<SimulationInfo> simInfo) = 0;
	virtual void makeFinalTestOutputHead(int numberOfPassedTests, int numberOfTests) = 0;
	virtual void makeFinalTestOutputFoot(int numberOfPassedTests, int numberOfTests) = 0;

private:

};
#endif