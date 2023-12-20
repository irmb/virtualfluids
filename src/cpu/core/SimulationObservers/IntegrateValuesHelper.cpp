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
#include "IntegrateValuesHelper.h"

#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbCuboid3D.h>
#include <vector>

#include "BCArray3D.h"
#include "BCSet.h"
#include "DataSet3D.h"
#include "LBMKernel.h"

//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::IntegrateValuesHelper(SPtr<Grid3D> grid, std::shared_ptr<vf::parallel::Communicator> comm, real minX1, real minX2,
                                             real minX3, real maxX1, real maxX2, real maxX3)
    :

      grid(grid), comm(comm), sVx1(0.0), sVx2(0.0), sVx3(0.0), sRho(0.0), sCellVolume(0.0), numberOfFluidsNodes(0),
      numberOfSolidNodes(0)
{
    boundingBox = GbCuboid3DPtr(new GbCuboid3D(minX1, minX2, minX3, maxX1, maxX2, maxX3));
    init(-1);
}
//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::IntegrateValuesHelper(SPtr<Grid3D> grid, std::shared_ptr<vf::parallel::Communicator> comm, real minX1, real minX2,
                                             real minX3, real maxX1, real maxX2, real maxX3, int level)
    :

      grid(grid), comm(comm), sVx1(0.0), sVx2(0.0), sVx3(0.0), sRho(0.0), sCellVolume(0.0), numberOfFluidsNodes(0),
      numberOfSolidNodes(0)
{
    boundingBox = GbCuboid3DPtr(new GbCuboid3D(minX1, minX2, minX3, maxX1, maxX2, maxX3));
    init(level);
}
//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::~IntegrateValuesHelper() = default;
//////////////////////////////////////////////////////////////////////////
void IntegrateValuesHelper::init(int level)
{
    using namespace vf::basics::constant;

    root = comm->isRoot();

    real orgX1, orgX2, orgX3;
    int gridRank = grid->getRank();
    int minInitLevel, maxInitLevel;
    if (level < 0) {
        minInitLevel = this->grid->getCoarsestInitializedLevel();
        maxInitLevel = this->grid->getFinestInitializedLevel();
    } else {
        minInitLevel = level;
        maxInitLevel = level;
    }

    real numSolids = c0o1;
    real numFluids = c0o1;
    for (int level_it = minInitLevel; level_it <= maxInitLevel; level_it++) {
        std::vector<SPtr<Block3D>> blockVector;
        grid->getBlocks(level_it, gridRank, blockVector);
        for (SPtr<Block3D> block : blockVector) {
            CalcNodes cn;
            cn.block = block;
            // Koords bestimmen
            UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

            orgX1 = val<1>(org);
            orgX2 = val<2>(org);
            orgX3 = val<3>(org);

            SPtr<LBMKernel> kernel                  = dynamicPointerCast<LBMKernel>(block->getKernel());
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            int ghostLayerWitdh                     = kernel->getGhostLayerWidth();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            real internX1, internX2, internX3;

            real dx               = grid->getDeltaX(block);
            UbTupleDouble3 orgDelta = grid->getNodeOffset(block);

            for (int ix3 = ghostLayerWitdh; ix3 < (int)distributions->getNX3() - ghostLayerWitdh; ix3++) {
                for (int ix2 = ghostLayerWitdh; ix2 < (int)distributions->getNX2() - ghostLayerWitdh; ix2++) {
                    for (int ix1 = ghostLayerWitdh; ix1 < (int)distributions->getNX1() - ghostLayerWitdh; ix1++) {
                        internX1 = orgX1 - val<1>(orgDelta) + ix1 * dx;
                        internX2 = orgX2 - val<2>(orgDelta) + ix2 * dx;
                        internX3 = orgX3 - val<3>(orgDelta) + ix3 * dx;
                        if (boundingBox->isPointInGbObject3D(internX1, internX2, internX3)) {
                            if (!bcArray->isSolid(ix1, ix2, ix3) && !bcArray->isUndefined(ix1, ix2, ix3)) {
                                cn.nodes.emplace_back(ix1, ix2, ix3);
                                numFluids++;
                            } else if (bcArray->isSolid(ix1, ix2, ix3)) {
                                numSolids++;
                            }
                        }
                    }
                }
            }
            if (cn.nodes.size() > 0)
                cnodes.push_back(cn);
        }
    }
    std::vector<real> rvalues;
    std::vector<real> values;
    values.push_back(numSolids);
    values.push_back(numFluids);
    rvalues = comm->gather(values);

    if (root) {
        numberOfSolidNodes  = c0o1;
        numberOfFluidsNodes = c0o1;
        int rsize           = (int)rvalues.size();
        int vsize           = (int)values.size();
        for (int i = 0; i < rsize; i += vsize) {
            numberOfSolidNodes += rvalues[i];
            numberOfFluidsNodes += rvalues[i + 1];
        }
    }
}
//////////////////////////////////////////////////////////////////////////
// calculation conventional rho, velocity and averaged data
void IntegrateValuesHelper::calculateAV()
{
    clearData();

    for (CalcNodes cn : cnodes) {
        SPtr<ILBMKernel> kernel                   = cn.block->getKernel();
        SPtr<AverageValuesArray3D> averagedValues = kernel->getDataSet()->getAverageValues();

        for (UbTupleInt3 node : cn.nodes) {
            real Avx = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVx);
            real Avy = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVy);
            real Avz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVz);

            real Avxx = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVxx);
            real Avyy = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVyy);
            real Avzz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVzz);

            real Avxz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVxz);
            sAvVx1 += std::abs(Avx);
            sAvVx2 += std::abs(Avy);
            sAvVx3 += std::abs(Avz);

            sTSx1 += sqrt(Avxx);
            sTSx2 += sqrt(Avyy);
            sTSx3 += sqrt(Avzz);

            sTSx1x3 += Avxz;
            numberOfFluidsNodes++;
        }
    }
    std::vector<real> values;
    std::vector<real> rvalues;
    values.push_back(sAvVx1);
    values.push_back(sAvVx2);
    values.push_back(sAvVx3);
    values.push_back(sTSx1);
    values.push_back(sTSx2);
    values.push_back(sTSx3);
    values.push_back(sTSx1x3);
    values.push_back(numberOfFluidsNodes);

    rvalues = comm->gather(values);
    if (root) {
        clearData();
        for (int i = 0; i < (int)rvalues.size(); i += 8) {
            sAvVx1 += rvalues[i];
            sAvVx2 += rvalues[i + 1];
            sAvVx3 += rvalues[i + 2];
            sTSx1 += rvalues[i + 3];
            sTSx2 += rvalues[i + 4];
            sTSx3 += rvalues[i + 5];
            sTSx1x3 += rvalues[i + 6];
            numberOfFluidsNodes += rvalues[i + 7];
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void IntegrateValuesHelper::calculateMQ()
{
    using namespace vf::basics::constant;

    real f[D3Q27System::ENDF + 1];
    real vx1, vx2, vx3, rho;
    clearData();

    // Funktionszeiger
    typedef void (*CalcMacrosFct)(const real *const & /*feq[27]*/, real & /*(d)rho*/, real & /*vx1*/,
                                  real & /*vx2*/, real & /*vx3*/);

    CalcMacrosFct calcMacros = NULL;

    for (CalcNodes cn : cnodes) {
        SPtr<ILBMKernel> kernel = cn.block->getKernel();
        real dx              = c1o1 / (real)(1 << cn.block->getLevel());
        real cellVolume      = dx * dx * dx;

        if (kernel->getCompressible()) {
            calcMacros = &D3Q27System::calcCompMacroscopicValues;
        } else {
            calcMacros = &D3Q27System::calcIncompMacroscopicValues;
        }

        SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
        SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
        for (UbTupleInt3 node : cn.nodes) {
            distributions->getPreCollisionDistribution(f, val<1>(node), val<2>(node), val<3>(node));
            calcMacros(f, rho, vx1, vx2, vx3);
            sRho += rho * cellVolume;
            sVx1 += vx1 * cellVolume;
            sVx2 += vx2 * cellVolume;
            sVx3 += vx3 * cellVolume;
            sCellVolume += cellVolume;
        }
    }
    std::vector<real> values(5);
    std::vector<real> rvalues;
    values[0] = sRho;
    values[1] = sVx1;
    values[2] = sVx2;
    values[3] = sVx3;
    values[4] = sCellVolume;

    rvalues = comm->gather(values);
    if (root) {
        clearData();
        int rsize = (int)rvalues.size();
        int vsize = (int)values.size();
        for (int i = 0; i < rsize; i += vsize) {
            sRho += rvalues[i];
            sVx1 += rvalues[i + 1];
            sVx2 += rvalues[i + 2];
            sVx3 += rvalues[i + 3];
            sCellVolume += rvalues[i + 4];
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void IntegrateValuesHelper::clearData()
{
    using namespace vf::basics::constant;
 
    sRho        = c0o1;
    sVx1        = c0o1;
    sVx2        = c0o1;
    sVx3        = c0o1;
    sCellVolume = c0o1;
    // sVm = 0.0;
    // sPress = 0.0;
    // numberOfFluidsNodes = 0.0;
    sAvVx1  = c0o1;
    sAvVx2  = c0o1;
    sAvVx3  = c0o1;
    sTSx1   = c0o1;
    sTSx2   = c0o1;
    sTSx3   = c0o1;
    sTSx1x3 = c0o1;
}
//////////////////////////////////////////////////////////////////////////
real IntegrateValuesHelper::getNumberOfFluidsNodes() { return this->numberOfFluidsNodes; }
//////////////////////////////////////////////////////////////////////////
real IntegrateValuesHelper::getNumberOfSolidNodes() { return this->numberOfSolidNodes; }
//////////////////////////////////////////////////////////////////////////
GbCuboid3DPtr IntegrateValuesHelper::getBoundingBox() { return this->boundingBox; }
//////////////////////////////////////////////////////////////////////////
std::vector<IntegrateValuesHelper::CalcNodes> IntegrateValuesHelper::getCNodes() { return cnodes; }

//! \}
