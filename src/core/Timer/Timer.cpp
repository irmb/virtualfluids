#include "Timer.h"
#include <ctime>


Timer::Timer()
{
    
}

Timer Timer::begin()
{
    Timer t;
    t.beginInClocks = clock();
    return t;
}

void Timer::end()
{
    const clock_t endInClocks = clock();
    timeInClocks = endInClocks - beginInClocks;
}

real Timer::getTimeInSeconds() const
{
    return real(timeInClocks) / CLOCKS_PER_SEC;
}