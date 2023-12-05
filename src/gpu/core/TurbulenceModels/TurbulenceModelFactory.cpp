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
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TurbulentViscosityFactory.cpp
//! \ingroup TurbulentViscosity
//! \author Henrik Asmuth
//=======================================================================================

#include "TurbulenceModelFactory.h"

#include <variant>

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include "PostProcessor/TurbulentViscosityKernels.h"

void TurbulenceModelFactory::setTurbulenceModel(vf::lbm::TurbulenceModel _turbulenceModel)
{
    this->turbulenceModel = _turbulenceModel;
    para->setTurbulenceModel(_turbulenceModel);
    if(this->turbulenceModel != vf::lbm::TurbulenceModel::None) para->setUseTurbulentViscosity(true);

    switch (this->turbulenceModel) {
        case vf::lbm::TurbulenceModel::AMD:
            this->turbulenceModelKernel = calcTurbulentViscosityAMD;
            break;
        default:
            this->turbulenceModelKernel = nullptr;
    }
}

void TurbulenceModelFactory::setModelConstant(real modelConstant)
{
    para->setSGSConstant(modelConstant);
}

void TurbulenceModelFactory::readConfigFile(const vf::basics::ConfigurationFile &configData)
{
    if (configData.contains("TurbulenceModel"))
    {
        std::string config = configData.getValue<std::string>("TurbulenceModel");
        
        if      (config == "Smagorinsky") this->setTurbulenceModel( vf::lbm::TurbulenceModel::Smagorinsky ); 
        else if (config == "AMD")         this->setTurbulenceModel( vf::lbm::TurbulenceModel::AMD );               
        else if (config == "QR" )         this->setTurbulenceModel( vf::lbm::TurbulenceModel::QR );             
        else if (config == "None")        this->setTurbulenceModel( vf::lbm::TurbulenceModel::None );           
        else    std::runtime_error("TurbulenceModelFactory: Invalid turbulence model!");           

        VF_LOG_INFO("Turbulence model: {}", config);
        
    }

    if (configData.contains("SGSconstant"))
    {
        para->setSGSConstant(configData.getValue<real>("SGSconstant"));

        VF_LOG_INFO("SGS constant: {}", para->getSGSConstant() );
    }
}

void TurbulenceModelFactory::runTurbulenceModelKernel(const int level) const
{
    if(this->turbulenceModelKernel) this->turbulenceModelKernel(para.get(), level);
}
