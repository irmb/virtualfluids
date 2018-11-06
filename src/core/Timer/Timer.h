#ifndef TIMER_H
#define TIMER_H

#include <chrono>

#include "VirtualFluidsDefinitions.h"

#include "DataTypes.h"

class VF_PUBLIC Timer
{
public:

    typedef std::chrono::high_resolution_clock::time_point timePoint;

    Timer();
    static Timer makeStart();

    void start();
    void end();

    real getTimeInSeconds() const;
    real getCurrentRuntimeInSeconds() const;

private:
    timePoint startTime;
    timePoint endTime;
};

#endif

