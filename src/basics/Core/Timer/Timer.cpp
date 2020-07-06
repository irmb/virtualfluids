#include "Timer.h"
#include "TimerImp.h"

#include <memory>

SPtr<Timer> Timer::makeStart()
{
    SPtr<Timer> t = std::make_shared<TimerImp>();
    t->start();
    return t;
}
