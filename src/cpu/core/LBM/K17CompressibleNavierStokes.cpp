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
//! \file K17CompressibleNavierStokes.cpp
//! \ingroup LBM
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================
#include "K17CompressibleNavierStokes.h"

#include <lbm/collision/CollisionParameter.h>
#include <lbm/collision/K17CompressibleNavierStokes.h>
#include <lbm/collision/TurbulentViscosity.h>

#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "EsoSplit.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "LBMKernel.h"

K17CompressibleNavierStokes::K17CompressibleNavierStokes()
{
    this->compressible = true;
}

void K17CompressibleNavierStokes::initDataSet()
{
    SPtr<DistributionArray3D> d(new EsoSplit(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    dataSet->setFdistributions(d);
}

SPtr<LBMKernel> K17CompressibleNavierStokes::clone()
{
    SPtr<LBMKernel> kernel(new K17CompressibleNavierStokes());
    kernel->setNX(nx);
    std::dynamic_pointer_cast<K17CompressibleNavierStokes>(kernel)->initDataSet();
    kernel->setCollisionFactor(this->collFactor);
    kernel->setBCSet(bcSet->clone(kernel));
    kernel->setWithForcing(withForcing);
    kernel->setForcingX1(muForcingX1);
    kernel->setForcingX2(muForcingX2);
    kernel->setForcingX3(muForcingX3);
    kernel->setIndex(ix1, ix2, ix3);
    kernel->setDeltaT(deltaT);
    kernel->setBlock(block.lock());

    return kernel;
}

void K17CompressibleNavierStokes::calculate(int step)
{
    (void)step; // parameter unused

    mu::value_type muX1, muX2, muX3;
    mu::value_type muDeltaT;
    mu::value_type muNu;

    if (withForcing) {
        muForcingX1.DefineVar("x1", &muX1); muForcingX1.DefineVar("x2", &muX2); muForcingX1.DefineVar("x3", &muX3);
        muForcingX2.DefineVar("x1", &muX1); muForcingX2.DefineVar("x2", &muX2); muForcingX2.DefineVar("x3", &muX3);
        muForcingX3.DefineVar("x1", &muX1); muForcingX3.DefineVar("x2", &muX2); muForcingX3.DefineVar("x3", &muX3);

        muDeltaT = deltaT;

        muForcingX1.DefineVar("dt", &muDeltaT);
        muForcingX2.DefineVar("dt", &muDeltaT);
        muForcingX3.DefineVar("dt", &muDeltaT);

        muNu = (c1o1 / c3o1) * (c1o1 / collFactor - c1o1 / c2o1);

        muForcingX1.DefineVar("nu", &muNu);
        muForcingX2.DefineVar("nu", &muNu);
        muForcingX3.DefineVar("nu", &muNu);
    }

    auto localDistributions = std::dynamic_pointer_cast<EsoSplit>(dataSet->getFdistributions())->getLocalDistributions();
    auto nonLocalDistributions = std::dynamic_pointer_cast<EsoSplit>(dataSet->getFdistributions())->getNonLocalDistributions();
    auto restDistributions = std::dynamic_pointer_cast<EsoSplit>(dataSet->getFdistributions())->getZeroDistributions();

    SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

    const int bcArrayMaxX1 = (int)bcArray->getNX1();
    const int bcArrayMaxX2 = (int)bcArray->getNX2();
    const int bcArrayMaxX3 = (int)bcArray->getNX3();

    const int minX1 = ghostLayerWidth;
    const int minX2 = ghostLayerWidth;
    const int minX3 = ghostLayerWidth;
    const int maxX1 = bcArrayMaxX1 - ghostLayerWidth;
    const int maxX2 = bcArrayMaxX2 - ghostLayerWidth;
    const int maxX3 = bcArrayMaxX3 - ghostLayerWidth;

    const real omega = collFactor;

    for (int x3 = minX3; x3 < maxX3; x3++) {
        for (int x2 = minX2; x2 < maxX2; x2++) {
            for (int x1 = minX1; x1 < maxX1; x1++) {

                if (bcArray->isSolid(x1, x2, x3) || bcArray->isUndefined(x1, x2, x3)) {
                    continue;
                }

                vf::lbm::CollisionParameter parameter;
                dataSet->getFdistributions()->getPreCollisionDistribution(parameter.distribution, x1, x2, x3);

                real forces[3] = { c0o1, c0o1, c0o1 };
                if (withForcing) // TODO: add level factor?
                {
                    muX1 = static_cast<real>(x1 - 1 + ix1 * maxX1);
                    muX2 = static_cast<real>(x2 - 1 + ix2 * maxX2);
                    muX3 = static_cast<real>(x3 - 1 + ix3 * maxX3);

                    const real forcingX1 = muForcingX1.Eval();
                    const real forcingX2 = muForcingX2.Eval();
                    const real forcingX3 = muForcingX3.Eval();

                    forces[0] = forcingX1 * deltaT;
                    forces[1] = forcingX2 * deltaT;
                    forces[2] = forcingX3 * deltaT;
                }
                parameter.forceX = forces[0];
                parameter.forceY = forces[1];
                parameter.forceZ = forces[2];

                parameter.omega = omega;

                real quadricLimiter[3] = { 0.01, 0.01, 0.01 }; // TODO: Where do we configure the quadricLimiter?
                parameter.quadricLimiter = quadricLimiter;

                vf::lbm::MacroscopicValues mv;  // not used
                vf::lbm::TurbulentViscosity tv; // not used
                vf::lbm::runK17CompressibleNavierStokes<vf::lbm::TurbulenceModel::None>(parameter, mv, tv);
                dataSet->getFdistributions()->setPostCollisionDistribution(parameter.distribution, x1, x2, x3);
            }
        }
    }
}
