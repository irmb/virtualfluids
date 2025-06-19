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
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
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

    //! \brief calculate the state of the next time step of the advection diffusion distributions
    void runADcollisionKernel(int level) const;

private:
    SPtr<Parameter> para;
    std::vector<SPtr<AdvectionDiffusionKernel>> adkernels;
};

#endif

//! \}
