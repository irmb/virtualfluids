#ifndef COLLISONSTRATEGY_H
#define COLLISONSTRATEGY_H

#include "UpdateGrid27.h"

std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t)>
    getFunctionForCollisionAndExchange(const bool useStreams, const int numberOfMpiProcesses,
                                       const bool kernelNeedsFluidNodeIndicesToRun);

class CollisionAndExchange_noStreams_indexKernel
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t);
};

class CollisionAndExchange_noStreams_oldKernel
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t);
};

class CollisionAndExchange_streams
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t);
};

#endif
