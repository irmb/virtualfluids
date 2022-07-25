#ifndef REFINEMENTSTRATEGY_H
#define REFINEMENTSTRATEGY_H

#include "UpdateGrid27.h"

//! \brief get a function which performs the interpolation between grid levels and performs the communication between gpus/ processes
//! \return a function to perform the interpolation and for multi-gpu simulations also the communication
std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level)>
    getFunctionForRefinementAndExchange(const bool useStreams, const int numberOfMpiProcesses, const int maxLevel,
                                        const bool useReducedCommunicationAfterFtoC) noexcept;

//! \brief Version of refinement: for multi-gpu simulations, with communication hiding ("streams"), only exchange the interpolated cells
//! \details recommended for multi-gpu simulations if the chosen collision kernel supports the use of cuda streams
class RefinementAndExchange_streams_exchangeInterface
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

//! \brief Version of refinement: for multi-gpu simulations, with communication hiding ("streams"), exchange all nodes
class RefinementAndExchange_streams_exchangeAllNodes
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

//! \brief Version of refinement: for multi-gpu simulations, without communication hiding ("streams"), only exchange the interpolated cells
//! \details recommended for multi-gpu simulations if the chosen collision kernel does NOT support the use of cuda streams
class RefinementAndExchange_noStreams_exchangeInterface
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

//! \brief Version of refinement: for multi-gpu simulations, without communication hiding ("streams"), exchange all nodes
class RefinementAndExchange_noStreams_exchangeAllNodes
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

//! Version of refinement: for single-gpu simulations
class Refinement_noExchange
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

//! Version of refinement: for uniform simulations (no grid refinement)
class NoRefinement
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

#endif
