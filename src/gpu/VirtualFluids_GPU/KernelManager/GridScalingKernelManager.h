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
//! \file GridScalingKernelManager.h
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GridScalingKernelManager_H
#define GridScalingKernelManager_H

#include "LBM/LB.h"
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"
#include "logger/Logger.h"
#include <functional>
#include <memory>
#include <stdexcept>

class Parameter;
class CudaMemoryManager;
class GridScalingFactory;
enum class CudaStreamIndex;
struct LBMSimulationParameter;
struct CUstream_st;

using gridScalingFC =
    std::function<void(LBMSimulationParameter *, LBMSimulationParameter *, ICellFC *, ICellNeighFC &, CUstream_st *stream)>;
using gridScalingCF =
    std::function<void(LBMSimulationParameter *, LBMSimulationParameter *, ICellCF *, ICellNeighCF &, CUstream_st *stream)>;

//! \class GridScalingKernelManager
//! \brief manage the cuda kernel calls
class VIRTUALFLUIDS_GPU_EXPORT GridScalingKernelManager
{
public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    //! \throws std::runtime_error when the user forgets to specify a scaling function
    GridScalingKernelManager(SPtr<Parameter> parameter, GridScalingFactory *gridScalingFactory);

    //! \brief calls the device function of the fine to coarse grid interpolation kernelH
    void runFineToCoarseKernelLB(const int level, InterpolationCellFineToCoarse *icellFC, ICellNeighFC &offFC, CudaStreamIndex streamIndex) const;

    //! \brief calls the device function of the fine to coarse grid interpolation kernel (advection diffusion)
    void runFineToCoarseKernelAD(const int level) const;

    //! \brief calls the device function of the coarse to fine grid interpolation kernel
    void runCoarseToFineKernelLB(const int level, InterpolationCellCoarseToFine *icellCF, ICellNeighCF &offCF, CudaStreamIndex streamIndex) const;

    //! \brief calls the device function of the coarse to fine grid interpolation kernel (advection diffusion)
    void runCoarseToFineKernelAD(const int level) const;

private:
    //! \brief check if grid scaling was set
    //! \throws std::runtime_error if interpolation nodes were assigned, but no scaling function was set in the grid
    //! scaling factory \param scalingFunctionFC: a kernel function for the grid scaling \param scalingStruct: a struct
    //! containing the grid nodes which are part of the interpolation \param scalingName: the name of the checked
    //! scaling function
    void checkScalingFunction(const gridScalingFC &scalingFunctionFC, const InterpolationCellFineToCoarse &scalingStruct,
                              const std::string &scalingName)
    {
        if (!scalingFunctionFC && scalingStruct.kFC > 0)
            throw std::runtime_error("The scaling function " + scalingName + " was not set!");
        if (scalingFunctionFC && scalingStruct.kFC == 0)
            VF_LOG_WARNING("The scaling function {} was set, although there is no refinement", scalingName);
    }

    //! \brief check if grid scaling was set
    //! \throws std::runtime_error if interpolation nodes were assigned, but no scaling function was set in the grid
    //! scaling factory \param scalingFunctionCF: a kernel function for the grid scaling \param scalingStruct: a struct
    //! containing the grid nodes which are part of the interpolation \param scalingName: the name of the checked
    //! scaling function
    void checkScalingFunction(const gridScalingCF &scalingFunctionCF, const InterpolationCellCoarseToFine &scalingStruct,
                              const std::string &scalingName)
    {
        if (!scalingFunctionCF && scalingStruct.kCF > 0)
            throw std::runtime_error("The scaling function " + scalingName + " was not set!");
        if (scalingFunctionCF && scalingStruct.kCF == 0)
            VF_LOG_WARNING("The scaling function {} was set, although there is no refinement", scalingName);
    }

    SPtr<Parameter> para;

    gridScalingFC scalingFineToCoarse = nullptr;
    gridScalingCF scalingCoarseToFine = nullptr;
};
#endif
