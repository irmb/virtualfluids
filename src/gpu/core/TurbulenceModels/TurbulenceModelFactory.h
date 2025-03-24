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
//! \author Henrik Asmuth
//=======================================================================================
#ifndef TurbulenceModelFactory_H
#define TurbulenceModelFactory_H

#include <functional>
#include <memory>
#include <optional>
#include <string>

namespace vf::basics
{
class ConfigurationFile;
}

#include <lbm/advectionDiffusion/TurbulentDiffusivity.h>
#include <lbm/collision/TurbulentViscosity.h>

class Parameter;

class TurbulenceModelFactory
{
public:
    TurbulenceModelFactory(std::shared_ptr<Parameter> parameter) : para(std::move(parameter))
    {
    }

    void setTurbulenceModel(std::string turbulenceModel);
    void setTurbulenceModel(vf::lbm::TurbulenceModel turbulenceModel);

    void setAdvectionDiffusionTurbulenceModel(std::string turbulenceModel);
    void setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel turbulenceModel);

    void setModelConstant(real modelConstant);

    void readConfigFile(const vf::basics::ConfigurationFile& configData);
    std::optional<std::function<void(Parameter*, int)>> getTurbulenceModelKernel()
    {
        return turbulenceModelKernel;
    }
    std::optional<std::function<void(Parameter*, int)>> getTurbulenceModelADKernel()
    {
        return turbulenceModelADKernel;
    }

private:
    std::optional<std::function<void(Parameter*, int)>> turbulenceModelKernel = std::nullopt;
    std::optional<std::function<void(Parameter*, int)>> turbulenceModelADKernel = std::nullopt;
    std::shared_ptr<Parameter> para;
};

#endif

//! \}
