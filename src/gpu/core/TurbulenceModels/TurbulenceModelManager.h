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
//! \addtogroup gpu_TurbulenceModels TurbulenceModels
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb
//=======================================================================================
#ifndef TurbulenceModelManager_H
#define TurbulenceModelManager_H
#include <functional>
#include <optional>
#include <memory>
#include "TurbulenceModelFactory.h"

namespace vf::gpu {

class Parameter;
class TurbulenceModelFactory;

class TurbulenceModelManager
{
public:
    TurbulenceModelManager(std::shared_ptr<Parameter> para, const std::shared_ptr<TurbulenceModelFactory>& factory) : para(std::move(para))
    {
        turbulenceModelKernel = factory->getTurbulenceModelKernel();
        turbulenceModelADKernel = factory->getTurbulenceModelADKernel();

    };

    void runTurbulenceModelKernel(int level) const;
    void runTurbulenceModelADKernel(int level) const;

private:
    std::shared_ptr<Parameter> para;
    std::optional<std::function<void(Parameter*, int)>> turbulenceModelKernel;
    std::optional<std::function<void(Parameter*, int)>> turbulenceModelADKernel;
};

}

#endif
//! \}
