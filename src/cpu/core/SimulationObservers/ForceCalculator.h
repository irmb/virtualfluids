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
//! \file ForceCalculator.h
//! \ingroup SimulationObservers
//! \author Sören Peters
//=======================================================================================

#include "lbm/constants/D3Q27.h"
 
#ifndef ForceCalculator_H
#define ForceCalculator_H

#include <memory>
#include <vector>

#include "Vector3D.h"

class D3Q27Interactor;
namespace vf::parallel {class Communicator;}
class DistributionArray3D;
class BoundaryConditions;

class ForceCalculator
{
public:
    ForceCalculator(std::shared_ptr<vf::parallel::Communicator> comm);
    virtual ~ForceCalculator();

    void calculateForces(std::vector<std::shared_ptr<D3Q27Interactor>> interactors);
    Vector3D getForces(int x1, int x2, int x3, std::shared_ptr<DistributionArray3D> distributions,
                       std::shared_ptr<BoundaryConditions> bc,
                       const Vector3D &boundaryVelocity = Vector3D(0.0, 0.0, 0.0)) const;

    Vector3D getGlobalForces() const;

private:
    void gatherGlobalForces();

    std::shared_ptr<vf::parallel::Communicator> comm;

    real forceX1global;
    real forceX2global;
    real forceX3global;
};

#endif
