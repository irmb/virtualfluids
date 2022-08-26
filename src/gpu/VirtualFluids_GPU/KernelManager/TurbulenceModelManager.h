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
//! \file TurbulenceModelManager.h
//! \ingroup KernelManager
//! \author Henrik Asmuth
//=======================================================================================
#ifndef TurbulenceModelManager_H
#define TurbulenceModelManager_H

#include <functional>
#include <memory>
#include <string>

#include "LBM/LB.h"
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

class CudaMemoryManager;
class TurbulenceModelFactory;
class Parameter;
struct LBMSimulationParameter;

using TurbulenceModelKernel = std::function<void(Parameter *, int )>;

//! \class BCKernelManager
//! \brief manage the cuda kernel calls to boundary conditions
//! \details This class stores the boundary conditions and manages the calls to the boundary condition kernels.
class VIRTUALFLUIDS_GPU_EXPORT TurbulenceModelManager
{
public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    //! \throws std::runtime_error when the user forgets to specify a boundary condition
    TurbulenceModelManager(SPtr<Parameter> parameter, TurbulenceModelFactory *turbulenceModelFactory);

    //! \brief calls the device function of the velocity boundary condition (post-collision)
    void runTurbulenceModelKernel(const int level) const;

private:
    //! \brief check if a boundary condition was set
    //! \throws std::runtime_error if boundary nodes were assigned, but no boundary condition was set in the boundary condition factory
    //! \param boundaryCondition: a kernel function for the boundary condition
    //! \param bcStruct: a struct containing the grid nodes which are part of the boundary condition
    //! \param bcName: the name of the checked boundary condition
   
    // template <typename bcFunction>
    // void checkTurbulenceModel(const bcFunction &boundaryCondition, const QforBoundaryConditions &bcStruct, const std::string &bcName)
    // {
    //     if (!boundaryCondition && bcStruct.numberOfBCnodes > 0)
    //         throw std::runtime_error("The boundary condition " + bcName + " was not set!");
    // }

    SPtr<Parameter> para;

    TurbulenceModelKernel turbulenceModelKernel = nullptr;

};
#endif
