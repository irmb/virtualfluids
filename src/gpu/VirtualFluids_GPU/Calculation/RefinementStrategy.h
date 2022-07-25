#ifndef REFINEMENTSTRATEGY_H
#define REFINEMENTSTRATEGY_H

#include "UpdateGrid27.h"

std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level)>
    getFunctionForRefinementAndExchange(const bool useStreams, const int numberOfMpiProcesses, const int maxLevel,
                                        const bool useReducedCommunicationAfterFtoC);

class RefinementAndExchange_streams_exchangeInterface
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

class RefinementAndExchange_streams_exchangeAllNodes
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

class RefinementAndExchange_noStreams_exchangeInterface
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

class RefinementAndExchange_noStreams_exchangeAllNodes
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

class Refinement_noExchange
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

class NoRefinement
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

#endif
