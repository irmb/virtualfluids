#ifndef TEST_H
#define TEST_H

#include "SimulationObserver.h"
#include "TestStatus.h"

#include <vector>

class Test : public SimulationObserver 
{
public:
	virtual void update() = 0;
	virtual void setSimulationCrashed() = 0;

	virtual TestStatus getTestStatus() = 0;
	virtual void makeConsoleOutput() = 0;

};
#endif 