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
#include <string>

#include <basics/PointerDefinitions.h>

#include "Calculation/Calculation.h"

class CudaMemoryManager;
class BoundaryConditionFactory;
class Parameter;
struct LBMSimulationParameter;

using boundaryCondition = std::function<void(LBMSimulationParameter *, QforBoundaryConditions *)>;
using boundaryConditionWithParameter = std::function<void(Parameter *, QforBoundaryConditions *, const int level)>;
using precursorBoundaryCondition = std::function<void(LBMSimulationParameter *, QforPrecursorBoundaryConditions *, real tRatio, real velocityRatio)>;

//! \class BCKernelManager
//! \brief manage the cuda kernel calls to boundary conditions
//! \details This class stores the boundary conditions and manages the calls to the boundary condition kernels.
class BoundaryConditionKernelManager
{
public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    //! \throws std::runtime_error when the user forgets to specify a boundary condition
    BoundaryConditionKernelManager(SPtr<Parameter> parameter, BoundaryConditionFactory *bcFactory);

    //! \brief calls the device function of the velocity boundary condition (post-collision)
    void runVelocityBCKernelPost(const int level) const;

    //! \brief calls the device function of the velocity boundary condition (pre-collision)
    void runVelocityBCKernelPre(const int level) const;

    //! \brief calls the device function of the geometry boundary condition (post-collision)
    void runGeoBCKernelPost(const int level) const;

    //! \brief calls the device function of the geometry boundary condition (pre-collision)
    void runGeoBCKernelPre(const int level, unsigned int t, CudaMemoryManager *cudaMemoryManager) const;

    //! \brief calls the device function of the slip boundary condition (post-collision)
    void runSlipBCKernelPost(const int level) const;

    //! \brief calls the device function of the no-slip boundary condition (post-collision)
    void runNoSlipBCKernelPost(const int level) const;

    //! \brief calls the device function of the pressure boundary condition (pre-collision)
    void runPressureBCKernelPre(const int level) const;

    //! \brief calls the device function of the precursor boundary condition
    void runPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager);

    //! \brief calls the device function of the outflow boundary condition
    void runOutflowBCKernelPre(const int level) const;

    //! \brief calls the device function of the stress wall model (post-collision)
    void runStressWallModelKernelPost(const int level) const;

private:
    //! \brief check if a boundary condition was set
    //! \throws std::runtime_error if boundary nodes were assigned, but no boundary condition was set in the boundary condition factory
    //! \param boundaryCondition: a kernel function for the boundary condition
    //! \param bcStruct: a struct containing the grid nodes which are part of the boundary condition
    //! \param bcName: the name of the checked boundary condition
    template <typename bcFunction, typename QforBC>
    void checkBoundaryCondition(const bcFunction &boundaryCondition, const QforBC &bcStruct, const std::string &bcName)
    {
        if (!boundaryCondition && bcStruct.numberOfBCnodes > 0)
            throw std::runtime_error("The boundary condition " + bcName + " was not set!");
    }

    void runDistributionPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager);
    void runVelocityPrecursorBCKernelPost(int level, uint t, CudaMemoryManager* cudaMemoryManager);

    SPtr<Parameter> para;

    boundaryCondition velocityBoundaryConditionPost = nullptr;
    boundaryCondition noSlipBoundaryConditionPost = nullptr;
    boundaryCondition slipBoundaryConditionPost = nullptr;
    boundaryCondition pressureBoundaryConditionPre = nullptr;
    boundaryCondition geometryBoundaryConditionPost = nullptr;
    boundaryConditionWithParameter stressBoundaryConditionPost = nullptr;
    precursorBoundaryCondition precursorBoundaryConditionPost = nullptr;
};
#endif

//! \}
