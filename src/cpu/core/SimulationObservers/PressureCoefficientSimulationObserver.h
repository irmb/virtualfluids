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
//! \file PressureCoefficientSimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef PressureCoefficientSimulationObserver_h__
#define PressureCoefficientSimulationObserver_h__

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"

class GbCuboid3D;
class D3Q27Interactor;
namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;

class PressureCoefficientSimulationObserver : public SimulationObserver
{
public:
    PressureCoefficientSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<GbCuboid3D> plane,
                                   const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm);
    ~PressureCoefficientSimulationObserver() override;

    void update(real step) override;

    void addInteractor(SPtr<D3Q27Interactor> interactor);
    void readValues(int step);

protected:
    void collectData(real step);
    void calculateRho();
    void writeValues(int step);

private:
    SPtr<GbCuboid3D> plane;
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;
    std::vector<SPtr<D3Q27Interactor>> interactors;
    int numberOfSteps;
    real maxStep;

    std::vector<UbTupleFloat3> nodes;
    std::vector<std::string> datanames;
    std::vector<std::vector<double>> data;

    std::vector<real> outValues;

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
};

#endif // PressureDistributionSimulationObserver_h__
