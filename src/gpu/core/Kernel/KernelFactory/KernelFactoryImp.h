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
#ifndef KERNEL_FACTORY_IMP_H
#define KERNEL_FACTORY_IMP_H

#include "KernelFactory.h"

class PorousMedia;

class KernelFactoryImp : public KernelFactory
{
public:
    std::vector< std::shared_ptr< Kernel>> makeKernels(std::shared_ptr<Parameter> para);
    std::vector< std::shared_ptr< AdvectionDiffusionKernel>> makeAdvDifKernels(std::shared_ptr<Parameter> para);

    void setPorousMedia(std::vector<std::shared_ptr<PorousMedia>> pm);

    std::shared_ptr<Kernel> makeKernel(std::shared_ptr<Parameter> para, std::string kernel, int level);
    std::shared_ptr<AdvectionDiffusionKernel> makeAdvDifKernel(std::shared_ptr<Parameter> para, std::string kernel, int level);
    void setKernelAtLevel(std::vector< std::shared_ptr<Kernel>> kernels, std::shared_ptr<Parameter> para, std::string kernel, int level);

    std::vector<std::shared_ptr<PorousMedia>> pm;
};
#endif
//! \}
