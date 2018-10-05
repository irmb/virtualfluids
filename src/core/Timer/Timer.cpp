#include "Timer.h"
#include <ctime>


Timer::Timer()
{
    
}

Timer Timer::makeStart()
{
    Timer t;
    t.start();
    return t;
}

void Timer::start()
{
    this->startTime = std::chrono::high_resolution_clock::now();
}

void Timer::end()
{
    this->endTime = std::chrono::high_resolution_clock::now();
}

real Timer::getTimeInSeconds() const
{
    return std::chrono::duration_cast<std::chrono::microseconds>( endTime - startTime ).count() / 1000000.0;
}

real Timer::getCurrentRuntimeInSeconds() const
{
    return std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::high_resolution_clock::now() - startTime ).count() / 1000000.0;
}
