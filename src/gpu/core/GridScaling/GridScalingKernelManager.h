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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GridScalingKernelManager_H
#define GridScalingKernelManager_H

#include "Calculation/Calculation.h"
#include <basics/PointerDefinitions.h>

#include <logger/Logger.h>
#include <functional>
#include <memory>
#include <stdexcept>

class Parameter;
class CudaMemoryManager;
class GridScalingFactory;
enum class CudaStreamIndex;
struct LBMSimulationParameter;
struct CUstream_st;

using gridScaling =
    std::function<void(LBMSimulationParameter *, LBMSimulationParameter *, ICells *, ICellNeigh &, CUstream_st *stream)>;

//! \class GridScalingKernelManager
//! \brief manage the cuda kernel calls
class GridScalingKernelManager
{
public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    //! \throws std::runtime_error when the user forgets to specify a scaling function
    GridScalingKernelManager(SPtr<Parameter> parameter, GridScalingFactory *gridScalingFactory);

    //! \brief calls the device function of the fine to coarse grid interpolation kernelH
    void runFineToCoarseKernelLB(const int level, InterpolationCells *fineToCoarse, ICellNeigh &neighborFineToCoarse, CudaStreamIndex streamIndex) const;

    //! \brief calls the device function of the coarse to fine grid interpolation kernel
    void runCoarseToFineKernelLB(const int level, InterpolationCells *coarseToFine, ICellNeigh &neighborCoarseToFine, CudaStreamIndex streamIndex) const;

private:
    //! \brief check if grid scaling was set
    //! \throws std::runtime_error if interpolation nodes were assigned, but no scaling function was set in the grid
    //! scaling factory \param scalingFunction: a kernel function for the grid scaling \param scalingStruct: a struct
    //! containing the grid nodes which are part of the interpolation \param scalingName: the name of the checked
    //! scaling function
    void checkScalingFunction(const gridScaling &scalingFunction, const InterpolationCells &scalingStruct,
                              const std::string &scalingName)
    {
        if (!scalingFunction && scalingStruct.numberOfCells > 0)
            throw std::runtime_error("The scaling function " + scalingName + " was not set!");
        if (scalingFunction && scalingStruct.numberOfCells == 0)
            VF_LOG_WARNING("The scaling function {} was set, although there is no refinement", scalingName);
    }

    SPtr<Parameter> para;

    gridScaling scalingFineToCoarse = nullptr;
    gridScaling scalingCoarseToFine = nullptr;
};
#endif
