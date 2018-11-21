#ifndef TEST_H
#define TEST_H

#include "SimulationObserver.h"

#include <memory>
#include <vector>
#include <string>

class TestSimulation;
class SimulationInfo;

class Test : public SimulationObserver 
{
public:
	virtual void update() = 0;
	virtual void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo) = 0;
	virtual std::string getLogFileOutput() = 0;
	virtual std::vector< bool> getPassedTests() = 0;
	virtual void makeConsoleOutput() = 0;

private:

};
#endif 