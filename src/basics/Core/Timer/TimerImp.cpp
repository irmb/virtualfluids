#include "TimerImp.h"


void TimerImp::start()
{
    this->startTime = std::chrono::high_resolution_clock::now();
}

void TimerImp::end()
{
    this->endTime = std::chrono::high_resolution_clock::now();
}

real TimerImp::getTimeInSeconds() const
{
    return std::chrono::duration_cast<std::chrono::microseconds>( endTime - startTime ).count() / 1000000.0;
}

real TimerImp::getCurrentRuntimeInSeconds() const
{
    return std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::high_resolution_clock::now() - startTime ).count() / 1000000.0;
}
