#ifndef TIMER_IMP_H
#define TIMER_IMP_H

#include "Timer.h"

#include <chrono>


#include "DataTypes.h"


class BASICS_EXPORT TimerImp : public Timer
{
public:
    using timePoint = std::chrono::high_resolution_clock::time_point;

    void start() override;
    void end() override;

    real getTimeInSeconds() const override;
    real getCurrentRuntimeInSeconds() const override;

private:
    timePoint startTime;
    timePoint endTime;
};

#endif

