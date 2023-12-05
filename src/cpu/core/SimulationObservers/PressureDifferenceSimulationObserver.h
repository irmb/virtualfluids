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
//! \file PressureDifferenceSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef D3Q27PRESSUREDIFFERENCESimulationObserver_H
#define D3Q27PRESSUREDIFFERENCESimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"
#include "LBMSystem.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class IntegrateValuesHelper;

class PressureDifferenceSimulationObserver : public SimulationObserver
{
public:
    PressureDifferenceSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                  SPtr<IntegrateValuesHelper> h1, SPtr<IntegrateValuesHelper> h2, real rhoReal,
                                  real uReal, real uLB,
                                  /*const SPtr<LBMUnitConverter> conv,*/ std::shared_ptr<vf::parallel::Communicator> comm);
    ~PressureDifferenceSimulationObserver() override;

    void update(real step) override;

protected:
    SPtr<IntegrateValuesHelper> h1, h2;
    std::string path;
    SPtr<LBMUnitConverter> conv;
    void collectData(real step);
    std::shared_ptr<vf::parallel::Communicator> comm;
    real factor1; //= (1/3)*rhoReal*(uReal/uLB)^2 for calculation pReal = rhoLB * (1/3)*rhoReal*(uReal/uLB)^2,
                     //rhoReal and uReal in SI
    real factor2; //= rhoReal*(uReal/uLB)^2       for calculation pReal = press * rhoReal*(uReal/uLB)^2, rhoReal and
                     //uReal in SI
};

#endif /* D3Q27RHODIFFERENCESimulationObserver_H_ */
