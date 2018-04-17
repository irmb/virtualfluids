#ifndef TIMER_H
#define TIMER_H

#include <VirtualFluidsDefinitions.h>
#include <DataTypes.h>

class VF_PUBLIC Timer
{
public:
    Timer();
    static Timer begin();
    void end();
    real getTimeInSeconds() const;

private:
    long beginInClocks;
    long timeInClocks;
};

#endif

