#ifndef COLLISONSTRATEGY_H
#define COLLISONSTRATEGY_H

#include "UpdateGrid27.h"

//! \brief get a function which performs the collision operator and performs the communication between gpus/ processes
//! \return a function to perform the collision and for multi-gpu simulations also the communication
std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t)>
    getFunctionForCollisionAndExchange(const bool useStreams, const int numberOfMpiProcesses,
                                       const bool kernelNeedsFluidNodeIndicesToRun);

//! \brief Version of collision: for multi-gpu simulations, without communication hiding ("streams"), for newer kernels that use an array of fluid nodes to determine which nodes to update
class CollisionAndExchange_noStreams_indexKernel
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t);
};

//! \brief Version of collision: for multi-gpu simulations, without communication hiding ("streams"), for old kernels
//! \details the only options for old kernel
class CollisionAndExchange_noStreams_oldKernel
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t);
};

//! \brief Version of collision: for multi-gpu simulations, with communication hiding ("streams"), for newer kernels that use an array of fluid nodes to determine which nodes to update
//! \details recommended for multi-gpu simulations if the chosen collision kernel supports the use of cuda streams
class CollisionAndExchange_streams
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t);
};

#endif
