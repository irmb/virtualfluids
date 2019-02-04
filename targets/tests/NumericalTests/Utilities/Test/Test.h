#ifndef TEST_H
#define TEST_H

#include "SimulationObserver.h"

#include <vector>


class Test : public SimulationObserver 
{
public:
	virtual void update() = 0;
	virtual std::vector<bool> getPassedTests() = 0;
	virtual void makeConsoleOutput() = 0;

private:

};
#endif 