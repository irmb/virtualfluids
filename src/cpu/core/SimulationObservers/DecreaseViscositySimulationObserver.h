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
//! \file DecreaseViscositySimulationObserver.h
//! \ingroup SimulationObservers
//! \author Sonja Uphoff
//=======================================================================================
#ifndef DecreaseViscositySimulationObserver_H
#define DecreaseViscositySimulationObserver_H

#include <PointerDefinitions.h>

#include "SimulationObserver.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"

#include "muParser.h"

class UbScheduler;
class Grid3D;
namespace vf::parallel {class Communicator;}

//! \brief The class sets viscosity/collision factor according to a previously defined function in time.
//! \details initialization in test case (example):
//! \code{.cpp}
//! mu::Parser decrViscFunc;                       //define a mu-parser function
//! decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");   //this function is time-dependent, the viscosity decreases a 1/t^2
//! decrViscFunc.DefineConst("nue0", nueLB);
//! decrViscFunc.DefineConst("c0", 0.1);           //constants such as c0 controll how fast the viscosity decreasis
//! SPtr<UbScheduler> DecrViscSch(new UbScheduler()); //the SimulationObserver is called according to a Scheduler
//! DecrViscSch->addSchedule(10,10,1000);          //in this case the viscosity is reset every 10 timesteps for the
//! first 1000 timesteps DecreaseViscositySimulationObserver decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm); \endcode
//! \author Sonja Uphoff

class DecreaseViscositySimulationObserver : public SimulationObserver
{
public:
    DecreaseViscositySimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, mu::Parser *nueFunc, std::shared_ptr<vf::parallel::Communicator> comm);
    ~DecreaseViscositySimulationObserver() override;
    //! calls collect PostprocessData.
    void update(real step) override;

protected:
    //! resets the collision factor depending on the current timestep.
    void setViscosity(real step);
    std::shared_ptr<vf::parallel::Communicator> comm;

private:
    mutable mu::value_type timeStep;
    mu::Parser *nueFunc;
};

#endif /* DecreaseViscositySimulationObserver_H_ */
