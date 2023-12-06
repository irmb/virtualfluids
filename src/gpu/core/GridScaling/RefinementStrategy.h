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
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Anna Wellmann, Martin Schönherr
//=======================================================================================
#ifndef REFINEMENTSTRATEGY_H
#define REFINEMENTSTRATEGY_H

#include <functional>

class UpdateGrid27;
class Parameter;

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

//! \brief Version of refinement: for single-gpu simulations
class Refinement_noExchange
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

//! \brief Version of refinement: for uniform simulations (no grid refinement)
class NoRefinement
{
public:
    void operator()(UpdateGrid27 *updateGrid, Parameter *para, int level);
};

#endif
