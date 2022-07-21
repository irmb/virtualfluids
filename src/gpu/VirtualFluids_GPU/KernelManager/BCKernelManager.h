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
//! \file BCKernelManager.h
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#ifndef BCKernelManager_H
#define BCKernelManager_H

#include <functional>
#include <memory>
#include <string>

#include "LBM/LB.h"
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

class CudaMemoryManager;
class BoundaryConditionFactory;
class Parameter;
struct LBMSimulationParameter;

using boundaryCondition = std::function<void(LBMSimulationParameter *, QforBoundaryConditions *)>;
using boundaryConditionPara = std::function<void(Parameter *, QforBoundaryConditions *, const int level)>;

//! \class BCKernelManager
//! \brief manage the cuda kernel calls to boundary conditions
class VIRTUALFLUIDS_GPU_EXPORT BCKernelManager
{
public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    BCKernelManager(SPtr<Parameter> parameter, BoundaryConditionFactory *bcFactory);

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

    //! \brief calls the device function of the pressure boundary condition (post-collision)
    void runPressureBCKernelPost(const int level) const;

    //! \brief calls the device function of the outflow boundary condition (pre-collision)
    void runOutflowBCKernelPre(const int level) const;

    //! \brief calls the device function of the stress wall model (post-collision)
    void runStressWallModelKernelPost(const int level) const;

private:
    void checkBoundaryCondition(const boundaryCondition &bc, const QforBoundaryConditions &bcStruct,
                                const std::string &bcName);
    void checkBoundaryCondition(const boundaryConditionPara &bc, const QforBoundaryConditions &bcStruct,
                                const std::string &bcName);

    SPtr<Parameter> para;

    boundaryCondition velocityBoundaryConditionPost;
    boundaryCondition noSlipBoundaryConditionPost;
    boundaryCondition slipBoundaryConditionPost;
    boundaryCondition pressureBoundaryConditionPre;
    boundaryCondition geometryBoundaryConditionPost;
    boundaryConditionPara stressBoundaryConditionPost;
};
#endif
