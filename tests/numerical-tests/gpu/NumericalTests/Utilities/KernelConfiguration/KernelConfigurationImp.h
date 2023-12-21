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
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#ifndef KERNEL_CONFIGURATION_IMP_H
#define KERNEL_CONFIGURATION_IMP_H

#include "KernelConfiguration.h"

#include <memory>

class KernelConfigurationImp : public KernelConfiguration
{
public:
    std::string getMainKernel();
    bool getMultiKernelOn();
    std::vector<int> getMultiKernelLevel();
    std::vector<std::string> getMultiKernel();

    static std::shared_ptr<KernelConfigurationImp> getNewInstance(std::string kernel);
    static std::shared_ptr<KernelConfigurationImp> getNewInstance(std::string kernel, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernel);

private:
    KernelConfigurationImp(std::string kernel);
    KernelConfigurationImp(std::string kernel, std::vector<int> multiKernelLevel, std::vector<std::string> multiKernel);
    KernelConfigurationImp() {};

    std::string mainKernel;
    bool multiKernelOn;
    std::vector<int> multiKernelLevel;
    std::vector<std::string> multiKernel;
};
#endif 
//! \}
