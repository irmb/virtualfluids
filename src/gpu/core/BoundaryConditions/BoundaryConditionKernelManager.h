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
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================
#ifndef BoundaryConditionKernelManager_H
#define BoundaryConditionKernelManager_H

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

#include <basics/PointerDefinitions.h>

#include "Calculation/Calculation.h"

namespace vf::gpu {

class CudaMemoryManager;
class BoundaryConditionFactory;
class Parameter;
struct LBMSimulationParameter;

using BoundaryConditionKernel = std::function<void(LBMSimulationParameter*, QforBoundaryConditions*)>;
using DirectionalBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, QforDirectionalBoundaryCondition*)>;
using BoundaryConditionWithParameterKernel = std::function<void(Parameter*, QforBoundaryConditions*, const int level)>;
using PrecursorBoundaryConditionKernel =
    std::function<void(LBMSimulationParameter*, QforPrecursorBoundaryConditions*, real tRatio, real velocityRatio)>;
using ADNoFluxBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, AdvectionDiffusionNoFluxBoundaryConditions bcParams)>;
using ADFluxBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, AdvectionDiffusionFluxBoundaryConditions bcParams)>;
using ADDirichletBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, AdvectionDiffusionDirichletBoundaryConditions bcParams)>;
using ADNeumannBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, AdvectionDiffusionNeumannBoundaryConditions bcParams)>;
using DirectionalADBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, QforDirectionalADBoundaryCondition*)>;
//! \class BCKernelManager
//! \brief manage the cuda kernel calls to boundary conditions
//! \details This class stores the boundary conditions and manages the calls to the boundary condition kernels.
class BoundaryConditionKernelManager
{
public:
    //! Class constructor
    //! \param bcFactory access to boundary condition factory without transfer of ownership
    //! \throws std::runtime_error when the user forgets to specify a boundary condition
    BoundaryConditionKernelManager(SPtr<Parameter> parameter, const BoundaryConditionFactory* bcFactory);

    //! \brief calls the device function of the velocity boundary condition (post-collision)
    void runVelocityBCKernelPost(int level) const;

    //! \brief calls the device function of the velocity boundary condition (pre-collision)
    void runVelocityBCKernelPre(int level) const;

    //! \brief calls the device function of the geometry boundary condition (post-collision)
    void runGeoBCKernelPost(int level) const;

    //! \brief calls the device function of the geometry boundary condition (pre-collision)
    void runGeoBCKernelPre(int level, unsigned int t, CudaMemoryManager* cudaMemoryManager) const;

    //! \brief calls the device function of the slip boundary condition (post-collision)
    void runSlipBCKernelPost(int level) const;

    //! \brief calls the device function of the no-slip boundary condition (post-collision)
    void runNoSlipBCKernelPost(int level) const;

    //! \brief calls the device function of the pressure boundary condition (pre-collision)
    void runPressureBCKernelPre(int level) const;

    //! \brief calls the device function of the precursor boundary condition
    void runPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager);

    //! \brief calls the device function of the stress wall model (post-collision)
    void runStressWallModelKernelPost(int level) const;

    void runADNoFluxBCKernel(int level) const ;
    void runADFluxBCKernel(int level) const ;
    void runADDirichletBCKernel(int level) const ;
    void runADNeumannBCKernel(int level) const ;
    //! \brief calls the device function of the directional advection-diffusion (temperature) outflow BC (post-collision)
    void runADDirectionalBCKernel(int level) const ;
        //! \brief calls the device function of the surface layer boundary condition (post-collision)
    void runSurfaceLayerBCKernelPost(int level) const;
private:
    //! \brief check if a directional boundary condition was set
    //! \throws std::runtime_error if boundary nodes were assigned, but no boundary condition was set in the boundary condition factory
    //! \param boundaryCondition: a kernel function for the boundary condition
    //! \param bcVector: a vector containing the boundary condition structs
    //! \param bcName: the name of the checked boundary condition
    template <typename bcFunction>
    void checkBoundaryCondition(const bcFunction& boundaryCondition,
                                const std::vector<QforDirectionalBoundaryCondition>& bcVector, const std::string& bcName)
    {
        if (!boundaryCondition && !bcVector.empty())
            throw std::runtime_error("The boundary condition " + bcName + " was not set!");
    }

    //! \brief check if a directional AD boundary condition was set
    //! \throws std::runtime_error if boundary nodes were assigned, but no boundary condition was set in the boundary condition factory
    //! \param boundaryCondition: a kernel function for the boundary condition
    //! \param bcVector: a vector containing the directional AD boundary condition structs
    //! \param bcName: the name of the checked boundary condition
    template <typename bcFunction>
    void checkBoundaryCondition(const bcFunction& boundaryCondition,
                                const std::vector<QforDirectionalADBoundaryCondition>& bcVector, const std::string& bcName)
    {
        if (!boundaryCondition && !bcVector.empty())
            throw std::runtime_error("The boundary condition " + bcName + " was not set!");
    }

    //! \brief check if a boundary condition was set
    //! \throws std::runtime_error if boundary nodes were assigned, but no boundary condition was set in the boundary condition factory
    //! \param boundaryCondition: a kernel function for the boundary condition
    //! \param bcStruct: a struct containing the grid nodes which are part of the boundary condition
    //! \param bcName: the name of the checked boundary condition
    template <typename bcFunction, typename QforBC>
    void checkBoundaryCondition(const bcFunction& boundaryCondition, const QforBC& bcStruct, const std::string& bcName)
    {
        if (!boundaryCondition && bcStruct.numberOfBCnodes > 0)
            throw std::runtime_error("The boundary condition " + bcName + " was not set!");
    }

    void runDistributionPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager);
    void runVelocityPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager);

    SPtr<Parameter> para;

    BoundaryConditionKernel velocityBoundaryConditionPost = nullptr;
    BoundaryConditionKernel noSlipBoundaryConditionPost = nullptr;
    BoundaryConditionKernel slipBoundaryConditionPost = nullptr;
    BoundaryConditionKernel geometryBoundaryConditionPost = nullptr;
    BoundaryConditionKernel stressBoundaryConditionPost = nullptr;
    PrecursorBoundaryConditionKernel precursorBoundaryConditionPost = nullptr;
    BoundaryConditionKernel pressureBoundaryConditionPre = nullptr;
    DirectionalBoundaryConditionKernel directionalPressureBoundaryConditionPre = nullptr;
    ADNoFluxBoundaryConditionKernel ADNoFluxBoundaryConditionPost = nullptr;
    ADFluxBoundaryConditionKernel ADFluxBoundaryConditionPost = nullptr;
    ADDirichletBoundaryConditionKernel ADDirichletBoundaryConditionPost = nullptr;
    ADNeumannBoundaryConditionKernel ADNeumannBoundaryConditionPost = nullptr;
    DirectionalADBoundaryConditionKernel ADDirectionalBoundaryConditionPost = nullptr;
    BoundaryConditionKernel surfaceLayerBoundaryConditionPost = nullptr;

};

}

#endif

//! \}
