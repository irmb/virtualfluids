//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_Calculation Calculation
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
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

//! \}
