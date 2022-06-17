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
//! \file ADKernelManager.h
//! \ingroup ADKernelManager
//! \author Martin Schoenherr
//=======================================================================================
#ifndef ADVECTION_DIFFUSION_H
#define ADVECTION_DIFFUSION_H

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "VirtualFluids_GPU_export.h"

//! \brief Class forwarding for Parameter, CudaMemoryManager
class Parameter;
class CudaMemoryManager;

//! \class ADKernelManager
//! \brief manage the advection diffusion kernel calls
class VIRTUALFLUIDS_GPU_EXPORT ADKernelManager{

public:
    //! \brief makes an object of ADKernelManager
    //! \param para shared pointer to instance of class Parameter
    static SPtr<ADKernelManager> make(SPtr<Parameter> parameter);

    //! \brief initialize the advection diffusion distributions
    void initAD(int level);
    
    //! \brief set initial concentration values at all nodes
    //! \param cudaManager instance of class CudaMemoryManager
    void setInitialNodeValuesAD(int level, SPtr<CudaMemoryManager> cudaMemoryManager);
    
    //! \brief calculate the state of the next time step of the advection diffusion distributions
    void runADcollisionKernel(int level);

    //! \brief calls the device function of the geometry boundary condition for advection diffusion
    void runADgeometryBCKernel(int level);

    //! \brief calls the device function of the velocity boundary condition for advection diffusion
    void runADveloBCKernel(int level);

    //! \brief calls the device function of the slip boundary condition for advection diffusion
    void runADslipBCKernel(int level);

    //! \brief calls the device function of the pressure boundary condition for advection diffusion
    void runADpressureBCKernel(int level);
    
    //! \brief copy the concentration from device to host and writes VTK file with concentration
    //! \param cudaManager instance of class CudaMemoryManager
    void printAD(int level, SPtr<CudaMemoryManager> cudaMemoryManager);


private:
    //! Class constructor
    //! \param parameter shared pointer to instance of class Parameter
    ADKernelManager(SPtr<Parameter> parameter);
    //! Class copy constructor
    //! \param ADKernelManager is a reference to ADKernelManager object
    ADKernelManager(const ADKernelManager&);

    //! \property para is a shared pointer to an object of Parameter
    SPtr<Parameter> para;
};

#endif
