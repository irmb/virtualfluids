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
#ifndef ADVECTION_DIFFUSION_H
#define ADVECTION_DIFFUSION_H

#include <vector>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>


//! \brief Class forwarding for Parameter, CudaMemoryManager
class Parameter;
class CudaMemoryManager;
class AdvectionDiffusionKernel;

//! \class ADKernelManager
//! \brief manage the advection diffusion kernel calls
class ADKernelManager
{

public:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    ADKernelManager(SPtr<Parameter> parameter, std::vector<SPtr<AdvectionDiffusionKernel>>& adkernels);

    //! \brief set initial concentration values at all nodes
    //! \param cudaMemoryManager instance of class CudaMemoryManager
    void setInitialNodeValuesAD(const int level, SPtr<CudaMemoryManager> cudaMemoryManager) const;

    //! \brief calculate the state of the next time step of the advection diffusion distributions
    void runADcollisionKernel(const int level) const;

    //! \brief calls the device function of the geometry boundary condition for advection diffusion
    void runADgeometryBCKernel(const int level) const;

    //! \brief calls the device function of the velocity boundary condition for advection diffusion
    void runADDirichletBCKernel(const int level) const;

    //! \brief calls the device function of the slip boundary condition for advection diffusion
    void runADslipBCKernel(const int level) const;

    //! \brief copy the concentration from device to host and writes VTK file with concentration
    //! \param cudaMemoryManager instance of class CudaMemoryManager
    void printAD(const int level, SPtr<CudaMemoryManager> cudaMemoryManager) const;

private:
    SPtr<Parameter> para;
    std::vector<SPtr<AdvectionDiffusionKernel>> adkernels;
};

#endif
