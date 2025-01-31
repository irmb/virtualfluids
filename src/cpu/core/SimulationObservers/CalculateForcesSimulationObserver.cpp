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
#include "CalculateForcesSimulationObserver.h"
#include "BCSet.h"

#include "BCArray3D.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include <parallel/Communicator.h>
#include "D3Q27Interactor.h"
#include "DataSet3D.h"
#include "DistributionArray3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

CalculateForcesSimulationObserver::CalculateForcesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                       std::shared_ptr<vf::parallel::Communicator> comm, real v, real a)
    : SimulationObserver(grid, s), path(path), comm(comm), v(v), a(a), forceX1global(0), forceX2global(0), forceX3global(0)
{
    if (comm->getProcessID() == comm->getRoot()) {
        std::ofstream ostr;
        std::string fname = path;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string file_path = ub_system::getPathFromString(fname);
            if (file_path.size() > 0) {
                ub_system::makeDirectory(file_path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr << "step" << ";" << "Cx" << ";" << "Cy" << ";" << "Cz" << ";" << "Fx" << ";" << "Fy" << ";" << "Fz" << std::endl;
        ostr.close();
    }
}
//////////////////////////////////////////////////////////////////////////
CalculateForcesSimulationObserver::~CalculateForcesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void CalculateForcesSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "D3Q27ForcesSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesSimulationObserver::collectData(real step)
{
    calculateForces();

    if (comm->getProcessID() == comm->getRoot()) {
        int istep = static_cast<int>(step);
        std::ofstream ostr;
        std::string fname = path;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string path = ub_system::getPathFromString(fname);
            if (path.size() > 0) {
                ub_system::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        calculateCoefficients();
        ostr << istep << ";" << C1 << ";" << C2 << ";" << C3 << ";" << forceX1global << ";" << forceX2global << ";" << forceX3global;
        ostr << std::endl;
        ostr.close();
    }
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesSimulationObserver::calculateForces()
{
    using namespace  vf::basics::constant;

    forceX1global = c0o1;
    forceX2global = c0o1;
    forceX3global = c0o1;

    for (SPtr<D3Q27Interactor> interactor : interactors) {
        for (BcNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap()) {
            real forceX1 = c0o1;
            real forceX2 = c0o1;
            real forceX3 = c0o1;

            SPtr<Block3D> block                             = t.first;
            std::set<std::vector<int>> &transNodeIndicesSet = t.second;

            SPtr<LBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            distributions->swap();

            int ghostLayerWidth = kernel->getGhostLayerWidth();
            int minX1           = ghostLayerWidth;
            int maxX1           = (int)bcArray->getNX1() - 1 - ghostLayerWidth;
            int minX2           = ghostLayerWidth;
            int maxX2           = (int)bcArray->getNX2() - 1 - ghostLayerWidth;
            int minX3           = ghostLayerWidth;
            int maxX3           = (int)bcArray->getNX3() - 1 - ghostLayerWidth;

            for (std::vector<int> node : transNodeIndicesSet) {
                int x1 = node[0];
                int x2 = node[1];
                int x3 = node[2];

                // without ghost nodes
                if (x1 < minX1 || x1 > maxX1 || x2 < minX2 || x2 > maxX2 || x3 < minX3 || x3 > maxX3)
                    continue;

                if (bcArray->isFluid(
                        x1, x2,
                        x3)) // es kann sein, dass der node von einem anderen interactor z.B. als solid gemarkt wurde!!!
                {
                    SPtr<BoundaryConditions> bc = bcArray->getBC(x1, x2, x3);
                    UbTupleDouble3 forceVec     = getForces(x1, x2, x3, distributions, bc);
                    forceX1 += val<1>(forceVec);
                    forceX2 += val<2>(forceVec);
                    forceX3 += val<3>(forceVec);
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
    std::vector<real> values;
    std::vector<real> rvalues;
    values.push_back(forceX1global);
    values.push_back(forceX2global);
    values.push_back(forceX3global);

    rvalues = comm->gather(values);
    if (comm->getProcessID() == comm->getRoot()) {
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
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 CalculateForcesSimulationObserver::getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions,
                                                     SPtr<BoundaryConditions> bc)
{
    using namespace  vf::basics::constant;

    UbTupleDouble3 force(c0o1, c0o1, c0o1);

    if (bc) {
        // references to tuple "force"
        double &forceX1 = val<1>(force);
        double &forceX2 = val<2>(force);
        double &forceX3 = val<3>(force);
        double f, fnbr;

        dynamicPointerCast<EsoTwist3D>(distributions)->swap();

        for (int fdir = d3q27_system::FSTARTDIR; fdir <= d3q27_system::FENDDIR; fdir++) {
            if (bc->hasNoSlipBoundaryFlag(fdir)) {
                const int invDir = d3q27_system::INVDIR[fdir];
                f = dynamicPointerCast<EsoTwist3D>(distributions)->getPostCollisionDistributionForDirection(x1, x2, x3, invDir);
                fnbr =
                    dynamicPointerCast<EsoTwist3D>(distributions)
                        ->getPostCollisionDistributionForDirection(x1 + d3q27_system::DX1[invDir], x2 + d3q27_system::DX2[invDir],
                                                         x3 + d3q27_system::DX3[invDir], fdir);

                forceX1 += (f + fnbr) * d3q27_system::DX1[invDir];
                forceX2 += (f + fnbr) * d3q27_system::DX2[invDir];
                forceX3 += (f + fnbr) * d3q27_system::DX3[invDir];
            }
        }
        dynamicPointerCast<EsoTwist3D>(distributions)->swap();
    }

    return force;
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesSimulationObserver::calculateCoefficients()
{
    using namespace  vf::basics::constant;
    
    real F1 = forceX1global;
    real F2 = forceX2global;
    real F3 = forceX3global;

    // return 2*F/(rho*v*v*a);
    C1 = c2o1 * F1 / (v * v * a);
    C2 = c2o1 * F2 / (v * v * a);
    C3 = c2o1 * F3 / (v * v * a);
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesSimulationObserver::addInteractor(SPtr<D3Q27Interactor> interactor) 
{ 
    interactors.push_back(interactor); 
}


//! \}
