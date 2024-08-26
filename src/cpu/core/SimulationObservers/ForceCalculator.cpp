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
//! \author Sören Peters, Konstantin Kutscher, Martin Geier
//=======================================================================================
#include "ForceCalculator.h"
#include "BCSet.h"

#include "BCArray3D.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include <parallel/Communicator.h>
#include "D3Q27Interactor.h"
#include "DataSet3D.h"
#include "DistributionArray3D.h"
#include "LBMKernel.h"
#include "EsoTwist3D.h"
#include <array>

ForceCalculator::ForceCalculator(std::shared_ptr<vf::parallel::Communicator> comm)
    : comm(comm), forceX1global(0), forceX2global(0), forceX3global(0)
{
}

ForceCalculator::~ForceCalculator() = default;

//////////////////////////////////////////////////////////////////////////
std::array<real, 3> ForceCalculator::getForces(int x1, int x2, int x3, std::shared_ptr<DistributionArray3D> distributions, std::shared_ptr<BoundaryConditions> bc) const
{
    using namespace  vf::basics::constant;

    real forceX1 = c0o1;
    real forceX2 = c0o1;
    real forceX3 = c0o1;

    if (bc) {
        dynamicPointerCast<EsoTwist3D>(distributions)->swap();

        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
            if (bc->hasNoSlipBoundaryFlag(fdir) || bc->hasVelocityBoundaryFlag(fdir)) {
                const int invDir = D3Q27System::INVDIR[fdir];
                const real f = dynamicPointerCast<EsoTwist3D>(distributions)->getPostCollisionDistributionForDirection(x1, x2, x3, invDir);
                const real fnbr = dynamicPointerCast<EsoTwist3D>(distributions)->getPostCollisionDistributionForDirection(x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

                forceX1 += (f + fnbr) * (real)D3Q27System::DX1[invDir];
                forceX2 += (f + fnbr) * (real)D3Q27System::DX2[invDir];
                forceX3 += (f + fnbr) * (real)D3Q27System::DX3[invDir];
            }
        }

        dynamicPointerCast<EsoTwist3D>(distributions)->swap();
    }

    return std::array<real, 3> { forceX1, forceX2, forceX3 };
}
//////////////////////////////////////////////////////////////////////////
std::array<real, 3> ForceCalculator::getForcesForDEM(int x1, int x2, int x3, std::shared_ptr<DistributionArray3D> distributions, std::shared_ptr<BoundaryConditions> bc, const std::array<real, 3>& boundaryVelocity) const
{
    using namespace vf::basics::constant;
    
    real forceX1 = c0o1;
    real forceX2 = c0o1;
    real forceX3 = c0o1;
    if (bc) {
        for (int fdir = d3q27_system::FSTARTDIR; fdir <= d3q27_system::FENDDIR; fdir++) {
            if (bc->hasNoSlipBoundaryFlag(fdir) || bc->hasVelocityBoundaryFlag(fdir)) {
                const int invDir  = d3q27_system::INVDIR[fdir];
                const real f    = distributions->getPostCollisionDistributionForDirection(x1, x2, x3, invDir);
                const real fnbr = distributions->getPostCollisionDistributionForDirection(
                    x1 + d3q27_system::DX1[invDir], x2 + d3q27_system::DX2[invDir], x3 + d3q27_system::DX3[invDir], fdir);

                real correction[3] = { c0o1, c0o1, c0o1 };
                if (bc->hasVelocityBoundaryFlag(fdir)) {
                    const real forceTerm = f - fnbr;
                    correction[0]          = forceTerm * boundaryVelocity[0];
                    correction[1]          = forceTerm * boundaryVelocity[1];
                    correction[2]          = forceTerm * boundaryVelocity[2];
                }

                // force consists of the MEM part and the Galilean invariance correction including the boundary velocity
                forceX1 += (f + fnbr) * d3q27_system::DX1[invDir] - correction[0];
                forceX2 += (f + fnbr) * d3q27_system::DX2[invDir] - correction[1];
                forceX3 += (f + fnbr) * d3q27_system::DX3[invDir] - correction[2];
            }
        }
    }
    return std::array<real, 3> { forceX1, forceX2, forceX3 };
}
//////////////////////////////////////////////////////////////////////////
std::array<real, 3> ForceCalculator::getForcesGalileanInvariantMomentumExchange(int x1, int x2, int x3, std::shared_ptr<DistributionArray3D> distributions, std::shared_ptr<BoundaryConditions> bc) const
{
    using namespace vf::basics::constant;

    real forceX1 = c0o1;
    real forceX2 = c0o1;
    real forceX3 = c0o1;

    if (bc) {
         dynamicPointerCast<EsoTwist3D>(distributions)->swap();

        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
            const int invDir = D3Q27System::INVDIR[fdir];
            if (bc->hasNoSlipBoundaryFlag(fdir) || bc->hasVelocityBoundaryFlag(fdir)) {
                real feq = c0o1, feqnbr = c0o1;
                if (bc->hasVelocityBoundaryFlag(fdir)) {
                    real vx1 = bc->getBoundaryVelocityX1();
                    real vx2 = bc->getBoundaryVelocityX2();
                    real vx3 = bc->getBoundaryVelocityX3();
                    feq = D3Q27System::getIncompFeqForDirection(fdir, 0.0, vx1, vx2, vx3);
                    feqnbr = D3Q27System::getIncompFeqForDirection(invDir, 0.0, vx1, vx2, vx3);
                }

                const real f = dynamicPointerCast<EsoTwist3D>(distributions)->getPostCollisionDistributionForDirection(x1, x2, x3, invDir);
                const real fnbr = dynamicPointerCast<EsoTwist3D>(distributions)->getPostCollisionDistributionForDirection(x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

                forceX1 += (f + fnbr - feq - feqnbr) * (real)D3Q27System::DX1[invDir];
                forceX2 += (f + fnbr - feq - feqnbr) * (real)D3Q27System::DX2[invDir];
                forceX3 += (f + fnbr - feq - feqnbr) * (real)D3Q27System::DX3[invDir];
            }
        }

        dynamicPointerCast<EsoTwist3D>(distributions)->swap();
    }

    return std::array<real, 3> { forceX1, forceX2, forceX3 };
}
//////////////////////////////////////////////////////////////////////////
void ForceCalculator::calculateForces(std::vector<SPtr<D3Q27Interactor>> interactors)
{
    using namespace vf::basics::constant;

    forceX1global = c0o1;
    forceX2global = c0o1;
    forceX3global = c0o1;

    for (const auto &interactor : interactors) {
        for (const auto &t : interactor->getBcNodeIndicesMap()) {
            real forceX1 = c0o1;
            real forceX2 = c0o1;
            real forceX3 = c0o1;

            SPtr<Block3D> block                     = t.first;
            SPtr<LBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            distributions->swap();

            const std::set<std::vector<int>> &transNodeIndices = t.second;
            for (std::vector<int> node : transNodeIndices) {
                int x1 = node[0];
                int x2 = node[1];
                int x3 = node[2];

                if (kernel->isInsideOfDomain(x1, x2, x3) && bcArray->isFluid(x1, x2, x3)) {
                    SPtr<BoundaryConditions> bc = bcArray->getBC(x1, x2, x3);
                    std::array<real, 3> forceVec = getForces(x1, x2, x3, distributions, bc);
                    forceX1 += forceVec[0];
                    forceX2 += forceVec[1];
                    forceX3 += forceVec[2];
                }
            }
            // if we have got discretization with more level
            // deltaX is LBM deltaX and equal LBM deltaT
            real deltaX = lbm_system::getDeltaT(block->getLevel()); // grid->getDeltaT(block);
            real deltaXquadrat = deltaX * deltaX;
            forceX1 *= deltaXquadrat;
            forceX2 *= deltaXquadrat;
            forceX3 *= deltaXquadrat;

            distributions->swap();

            forceX1global += forceX1;
            forceX2global += forceX2;
            forceX3global += forceX3;
        }
    }
    gatherGlobalForces();
}

void ForceCalculator::gatherGlobalForces()
{
    using namespace vf::basics::constant;

    std::vector<real>
        values; // intel compiler 17 dasn't support this { forceX1global , forceX2global, forceX3global };
    values.push_back(forceX1global);
    values.push_back(forceX2global);
    values.push_back(forceX3global);
    std::vector<real> rvalues = comm->gather(values);

    if (comm->isRoot()) {
        forceX1global = c0o1;
        forceX2global = c0o1;
        forceX3global = c0o1;

        for (int i = 0; i < (int)rvalues.size(); i += 3) {
            forceX1global += rvalues[i];
            forceX2global += rvalues[i + 1];
            forceX3global += rvalues[i + 2];
        }
    }
}

std::array<real, 3> ForceCalculator::getGlobalForces() const
{
    return std::array<real, 3> { forceX1global, forceX2global, forceX3global };
}

//! \}
