#ifndef TEST_H
#define TEST_H

#include "SimulationObserver.h"
#include "TestStatus.h"

#include <vector>
#include <string>

class Test : public SimulationObserver 
{
public:
    virtual ~Test() = default;
    virtual void run() = 0;

    virtual TestStatus getTestStatus() = 0;
    virtual void makeConsoleOutput() = 0;

};
#endif 