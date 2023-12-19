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
//! \addtogroup cpu_SimulationObservers SimulationObservers
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef D3Q27ForcesSimulationObserver_H
#define D3Q27ForcesSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "UbTuple.h"
#include "lbm/constants/D3Q27.h"

class ForceCalculator;
namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class DistributionArray3D;
class BoundaryConditions;

class CalculateForcesSimulationObserver : public SimulationObserver
{
public:
    //! Constructor
    //! \param v - velocity of fluid in LB units
    //! \param a - area of object in LB units
    CalculateForcesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm,
                               real v, real a);
    ~CalculateForcesSimulationObserver() override;
    void update(real step) override;
    void addInteractor(SPtr<D3Q27Interactor> interactor);

protected:
    void collectData(real step);
    void calculateForces();
    UbTupleDouble3 getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions,
                             SPtr<BoundaryConditions> bc);
    void calculateCoefficients();
    void write(std::ofstream *fileObject, real value, char *separator);

private:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;
    std::vector<SPtr<D3Q27Interactor>> interactors;
    real forceX1global;
    real forceX2global;
    real forceX3global;
    real v; //!< is the speed of the object relative to the fluid
    real a; //!< is the reference area
    real C1;
    real C2;
    real C3;
};

#endif /* D3Q27ForcesSimulationObserver_H */

//! \}
