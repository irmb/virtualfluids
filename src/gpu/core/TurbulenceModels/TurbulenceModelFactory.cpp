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

#include "TurbulenceModelFactory.h"

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

#include "Parameter/Parameter.h"
#include "PostProcessor/TurbulentViscosityKernels.h"
#include "lbm/collision/TurbulentViscosity.h"

using ADTurbulenceModel = vf::lbm::advection_diffusion::TurbulenceModel;


void TurbulenceModelFactory::setTurbulenceModel(std::string turbulenceModel)
{
    if (turbulenceModel == "Smagorinsky")
        this->setTurbulenceModel(vf::lbm::TurbulenceModel::Smagorinsky);
    else if (turbulenceModel == "AMD")
        this->setTurbulenceModel(vf::lbm::TurbulenceModel::AMD);
    else if (turbulenceModel == "AMDStratified")
        this->setTurbulenceModel(vf::lbm::TurbulenceModel::AMDStratified);
    else if (turbulenceModel == "QR")
        this->setTurbulenceModel(vf::lbm::TurbulenceModel::QR);
    else if (turbulenceModel == "None")
        this->setTurbulenceModel(vf::lbm::TurbulenceModel::None);
    else
        throw std::runtime_error("TurbulenceModelFactory: Invalid turbulence model! Model name found: " + turbulenceModel);

    VF_LOG_INFO("Turbulence model: {}", turbulenceModel);
}

void TurbulenceModelFactory::setTurbulenceModel(vf::lbm::TurbulenceModel turbulenceModel)
{
    para->setTurbulenceModel(turbulenceModel);

    if (turbulenceModel != vf::lbm::TurbulenceModel::None && para->getSGSConstant() == c0o1)
        throw std::runtime_error("Turbulence Model requires SGS constant!");

    if (turbulenceModel != vf::lbm::TurbulenceModel::None)
        para->setUseTurbulentViscosity(true);

    if (turbulenceModel == vf::lbm::TurbulenceModel::AMD)
        this->turbulenceModelKernel = calculateTurbulentViscosityAMD;

    if (turbulenceModel == vf::lbm::TurbulenceModel::AMDStratified)
        this->turbulenceModelKernel = calculateTurbulentViscosityAndDiffusivityAMDStratified;
}

void TurbulenceModelFactory::setModelConstant(real modelConstant)
{
    para->setSGSConstant(modelConstant);
}

void TurbulenceModelFactory::setAdvectionDiffusionTurbulenceModel(std::string turbulenceModel)
{
    VF_LOG_INFO("Turbulence Model Advection Diffuision {}", turbulenceModel);
    if (turbulenceModel == "Default") {
        this->setAdvectionDiffusionTurbulenceModel(ADTurbulenceModel::Default);
        VF_LOG_INFO("Turbulent Prandtl Number: {}", para->getTurbulentPrandtlNumber());
    } else if (turbulenceModel == "Moeng")
        this->setAdvectionDiffusionTurbulenceModel(ADTurbulenceModel::Moeng);
    else if (turbulenceModel == "AMDStratified")
        this->setAdvectionDiffusionTurbulenceModel(ADTurbulenceModel::AMDStratified);
    else if (turbulenceModel == "None")
        this->setAdvectionDiffusionTurbulenceModel(ADTurbulenceModel::None);
    else
        throw std::runtime_error("TurbulenceModelFactory: Invalid advection diffusion turbulence model!");
}

void TurbulenceModelFactory::setAdvectionDiffusionTurbulenceModel(ADTurbulenceModel turbulenceModel)
{
    if (!para->getUseTurbulentViscosity() && turbulenceModel != ADTurbulenceModel::None)
        throw std::runtime_error(
            "TurbulenceModelFactory: Turbulent viscosity must be enabled to use an advection diffusion turbulence model!");

    if (turbulenceModel == ADTurbulenceModel::Default && para->getTurbulentPrandtlNumber() == c0o1)
        throw std::runtime_error(
            "TurbulenceModelFactory: Prandtl number must be set to use the default advection diffusion turbulence model!");
        
    if((turbulenceModel == ADTurbulenceModel::AMDStratified) == (para->getTurbulenceModel() != vf::lbm::TurbulenceModel::AMDStratified))
        throw std::runtime_error("TurbulenceModelFactory: Can only use AMDstratified for turbulent viscosity and diffusivity together!");

    if(turbulenceModel != ADTurbulenceModel::None && turbulenceModel != ADTurbulenceModel::Default && para->getTurbulentPrandtlNumber() != c0o1)
        VF_LOG_INFO("Turbulent Prandtl Number is set but AD Turbulence Model does not use turbulent Prandtl Number");
    
    para->setAdvectionDiffusionTurbulenceModel(turbulenceModel);
    if (turbulenceModel != ADTurbulenceModel::None)
        para->setUseTurbulentDiffusivity(true);
    if (turbulenceModel == ADTurbulenceModel::Moeng)
        this->turbulenceModelADKernel = calculateTurbulentDiffusivityMoeng;
}

void TurbulenceModelFactory::readConfigFile(const vf::basics::ConfigurationFile& configData)
{
    const std::string SGSkey = "SGSconstant";
    const std::string turbulenceModelKey = "TurbulenceModel";
    const std::string ADTurbulenceModelKey = "TurbulenceModelAdvectionDiffusion";
    if (configData.contains(SGSkey)) {
        para->setSGSConstant(configData.getValue<real>(SGSkey));
        VF_LOG_INFO("SGS constant: {}", para->getSGSConstant());
    }

    if (configData.contains(turbulenceModelKey)) {
        const std::string config = configData.getValue<std::string>(turbulenceModelKey);
        setTurbulenceModel(config);
    }

    if (configData.contains(ADTurbulenceModelKey)) {
        const std::string config = configData.getValue<std::string>(ADTurbulenceModelKey);
        setAdvectionDiffusionTurbulenceModel(config);
    }
}

void TurbulenceModelFactory::runTurbulenceModelKernel(int level) const
{
    if (this->turbulenceModelKernel.has_value())
        this->turbulenceModelKernel.value()(para.get(), level);
}

void TurbulenceModelFactory::runTurbulenceModelADKernel(int level) const
{
    if (this->turbulenceModelADKernel.has_value())
        this->turbulenceModelADKernel.value()(para.get(), level);
}

//! \}
