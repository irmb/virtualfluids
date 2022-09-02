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
//! \file TurbulentViscosityFactory.h
//! \ingroup TurbulentViscosity
//! \author Henrik Asmuth
//=======================================================================================
#ifndef TurbulenceModelFactory_H
#define TurbulenceModelFactory_H

#include <functional>
#include <map>
#include <string>
#include <variant>

#include "LBM/LB.h"
#include "Parameter/Parameter.h"

#include <basics/config/ConfigurationFile.h>

class Parameter;

using TurbulenceModelKernel = std::function<void(Parameter *, int )>;

class TurbulenceModelFactory
{
public:
    
    TurbulenceModelFactory(SPtr<Parameter> parameter): para(parameter) {}

    void setTurbulenceModel(TurbulenceModel _turbulenceModel);

    void setModelConstant(real modelConstant);

    void readConfigFile(const vf::basics::ConfigurationFile &configData);

    void runTurbulenceModelKernel(const int level) const;

private:
    TurbulenceModel turbulenceModel = TurbulenceModel::None;
    TurbulenceModelKernel turbulenceModelKernel = nullptr;
    SPtr<Parameter> para;

};

#endif
