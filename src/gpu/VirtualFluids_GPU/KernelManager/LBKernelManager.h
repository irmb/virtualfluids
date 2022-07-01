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
//! \file LBKernelManager.h
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#ifndef LBKernelManager_H
#define LBKernelManager_H

#include <functional>
#include <memory>

#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"
#include "LBM/LB.h"

class CudaMemoryManager;
class BoundaryConditionFactory;
class Parameter;
struct LBMSimulationParameter;

using boundaryCondition = std::function<void(LBMSimulationParameter *, QforBoundaryConditions *)>;

//! \class LBKernelManager
//! \brief manage the cuda kernel calls
class VIRTUALFLUIDS_GPU_EXPORT LBKernelManager
{
public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    LBKernelManager(SPtr<Parameter> parameter, BoundaryConditionFactory *bcFactory);

    void setBoundaryConditionKernels();

    //! \brief calls the device function of the lattice Boltzmann kernel
    void runLBMKernel(const int level) const;

    //! \brief calls the device function of the velocity boundary condition (post-collision)
    void runVelocityBCKernelPost(const int level) const;

    //! \brief calls the device function of the velocity boundary condition (pre-collision)
    void runVelocityBCKernelPre(const int level) const;

    //! \brief calls the device function of the geometry boundary condition (post-collision)
    void runGeoBCKernelPost(const int level) const;

    //! \brief calls the device function of the geometry boundary condition (pre-collision)
    void runGeoBCKernelPre(const int level, unsigned int t, CudaMemoryManager *cudaMemoryManager) const;

    //! \brief calls the device function of the slip boundary condition
    void runSlipBCKernel(const int level) const;

    //! \brief calls the device function of the no-slip boundary condition
    void runNoSlipBCKernel(const int level) const;

    //! \brief calls the device function of the pressure boundary condition (pre-collision)
    void runPressureBCKernelPre(const int level) const;

    //! \brief calls the device function of the pressure boundary condition (post-collision)
    void runPressureBCKernelPost(const int level) const;

    //! \brief calls the device function of the outflow boundary condition
    void runOutflowBCKernelPre(const int level) const;

    //! \brief calls the device function of the stress wall model
    void runStressWallModelKernel(const int level) const;

    //! \brief calls the device function that calculates the macroscopic values
    void calculateMacroscopicValues(const int level) const;

private:
    SPtr<Parameter> para;

    boundaryCondition velocityBoundaryConditionPost;
};
#endif
