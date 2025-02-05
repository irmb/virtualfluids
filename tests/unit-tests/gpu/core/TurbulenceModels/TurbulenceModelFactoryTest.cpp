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
//! \addtogroup gpu_turbulenceModel_tests TurbulenceModelFactory
//! \ingroup gpu_core_tests core
//! \{
//! \author Henry Korb
//=======================================================================================
#include <gmock/gmock.h>

#include <basics/config/ConfigurationFile.h>
#include <gpu/core/TurbulenceModels/TurbulenceModelFactory.h>
#include <stdexcept>

#include "../testUtilitiesGPU.h"
#include "basics/tests/testUtilities.h"

#include "Parameter/Parameter.h"
#include "PostProcessor/TurbulentViscosityKernels.h"
#include "lbm/advectionDiffusion/TurbulentDiffusivity.h"
#include "lbm/collision/TurbulentViscosity.h"

class TurbulenceModelFactoryTest_Initialization : public testing::Test
{
protected:
    TurbulenceModelFactory* tmFactory = nullptr;
    SPtr<Parameter> para = std::make_shared<Parameter>();

    void SetUp() override
    {
        tmFactory = new TurbulenceModelFactory(para);
    }
};

using turbulenceModelFunction = void (*)(Parameter*, int);

turbulenceModelFunction getTurbulenceModelKernel(std::optional<std::function<void(Parameter*, int)>> kernelOption)
{
    if(!kernelOption)
        return nullptr;
    auto kernel = kernelOption.value();
    void (*tmTarget)(Parameter*, int) = (*kernel.target<void (*)(Parameter*, int)>());
    return tmTarget;
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_sgs)
{
    vf::basics::ConfigurationFile config;
    config.data["SGSconstant"] = "0.3333";
    tmFactory->readConfigFile(config);
    EXPECT_THAT(para->getSGSConstant(), RealEq(0.3333));
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_none)
{
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::None);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentViscosity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_smagorinsky_no_sgs)
{
    EXPECT_THROW(tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::Smagorinsky), std::runtime_error);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_smagorinsky_with_sgs)
{
    para->setSGSConstant(c1o1);
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::Smagorinsky);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::Smagorinsky);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_amd)
{
    para->setSGSConstant(c1o1);
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::AMD);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::AMD);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), calculateTurbulentViscosityAMD);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_qr)
{
    para->setSGSConstant(c1o1);
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::QR);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::QR);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_none_str)
{
    tmFactory->setTurbulenceModel("None");
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentViscosity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_smagorinsky_str)
{
    para->setSGSConstant(c1o1);
    tmFactory->setTurbulenceModel("Smagorinsky");
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::Smagorinsky);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_amd_str)
{
    para->setSGSConstant(c1o1);
    tmFactory->setTurbulenceModel("AMD");
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::AMD);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), calculateTurbulentViscosityAMD);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_qr_str)
{
    para->setSGSConstant(c1o1);
    tmFactory->setTurbulenceModel("QR");
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::QR);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_empty)
{
    vf::basics::ConfigurationFile config;
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentViscosity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_none)
{
    vf::basics::ConfigurationFile config;
    config.data["None"] = "None";
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentViscosity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_smagorinsky)
{
    vf::basics::ConfigurationFile config;
    config.data["TurbulenceModel"] = "Smagorinsky";
    config.data["SGSconstant"] = "0.33333";
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::Smagorinsky);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_amd)
{
    vf::basics::ConfigurationFile config;
    config.data["TurbulenceModel"] = "AMD";
    config.data["SGSconstant"] = "0.33333";

    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::AMD);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), calculateTurbulentViscosityAMD);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_qr)
{
    vf::basics::ConfigurationFile config;
    config.data["TurbulenceModel"] = "QR";
    config.data["SGSconstant"] = "0.33333";

    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::QR);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_none)
{
    tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::None);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}


TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_default_no_turb_visc)
{
    EXPECT_THROW(tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::Default), std::runtime_error);
}


TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_default_no_prandtl)
{
    para->setUseTurbulentViscosity(true);
    EXPECT_THROW(tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::Default), std::runtime_error);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_default)
{
    para->setUseTurbulentViscosity(true);
    para->setTurbulentPrandtlNumber(0.7);
    tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::Default);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::Default);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_moeng_no_turb_viscosity)
{
    EXPECT_THROW(tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::Moeng), std::runtime_error);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_moeng)
{
    para->setUseTurbulentViscosity(true);
    tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::Moeng);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::Moeng);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), calculateTurbulentDiffusivityMoeng);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AMD_Stratified_OnlyViscosity)
{
    para->setSGSConstant(0.3333);
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::AMDStratified);
    EXPECT_THROW(
    tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::None), std::runtime_error);
}
TEST_F(TurbulenceModelFactoryTest_Initialization, set_AMD_Stratified_OnlyDiffusivity)
{
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::None);
    para->setUseTurbulentViscosity(true);
    EXPECT_THROW(
    tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::AMDStratified), std::runtime_error);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AMD_Stratified)
{
    para->setSGSConstant(0.3333);
    tmFactory->setTurbulenceModel(vf::lbm::TurbulenceModel::AMDStratified);
    tmFactory->setAdvectionDiffusionTurbulenceModel(vf::lbm::advection_diffusion::TurbulenceModel::AMDStratified);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::AMDStratified);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), calculateTurbulentViscosityAndDiffusivityAMDStratified);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_none_str)
{
    tmFactory->setAdvectionDiffusionTurbulenceModel("None");
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_default_str)
{
    para->setUseTurbulentViscosity(true);
    para->setTurbulentPrandtlNumber(0.7);
    tmFactory->setAdvectionDiffusionTurbulenceModel("Default");
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::Default);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, set_AD_moeng_Str)
{
    para->setUseTurbulentViscosity(true);
    tmFactory->setAdvectionDiffusionTurbulenceModel("Moeng");
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::Moeng);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), calculateTurbulentDiffusivityMoeng);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_AD_none)
{
    vf::basics::ConfigurationFile config;
    config.data["TurbulenceModelAdvectionDiffusion"] = "None";
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::None);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), false);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_AD_default)
{
    vf::basics::ConfigurationFile config;
    config.data["TurbulenceModelAdvectionDiffusion"] = "Default";
    para->setUseTurbulentViscosity(true);
    para->setTurbulentPrandtlNumber(0.7);
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::Default);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_AD_moeng)
{
    vf::basics::ConfigurationFile config;
    config.data["TurbulenceModelAdvectionDiffusion"] = "Moeng";
    para->setUseTurbulentViscosity(true);
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::Moeng);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), calculateTurbulentDiffusivityMoeng);
}

TEST_F(TurbulenceModelFactoryTest_Initialization, read_amd_stratified)
{
    vf::basics::ConfigurationFile config;
    config.data["SGSconstant"] = "0.333";
    config.data["TurbulenceModel"] = "AMDStratified";
    config.data["TurbulenceModelAdvectionDiffusion"] = "AMDStratified";
    tmFactory->readConfigFile(config);
    EXPECT_EQ(para->getADTurbulenceModel(), vf::lbm::advection_diffusion::TurbulenceModel::AMDStratified);
    EXPECT_EQ(para->getTurbulenceModel(), vf::lbm::TurbulenceModel::AMDStratified);
    EXPECT_EQ(para->getUseTurbulentDiffusivity(), true);
    EXPECT_EQ(para->getUseTurbulentViscosity(), true);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelKernel()), calculateTurbulentViscosityAndDiffusivityAMDStratified);
    EXPECT_EQ(getTurbulenceModelKernel(tmFactory->getTurbulenceModelADKernel()), nullptr);
}