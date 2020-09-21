#ifndef TIMER_H
#define TIMER_H

#include "basics_export.h"

#include "DataTypes.h"
#include "PointerDefinitions.h"

class BASICS_EXPORT Timer
{
public:

    static SPtr<Timer> makeStart();

    virtual void start() = 0;
    virtual void end() = 0;

    virtual real getTimeInSeconds() const = 0;
    virtual real getCurrentRuntimeInSeconds() const = 0;
};

#endif

