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
#include "TurbulenceIntensitySimulationObserver.h"

#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "DataSet3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"
#include "basics/utilities/UbMath.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

TurbulenceIntensitySimulationObserver::TurbulenceIntensitySimulationObserver(SPtr<Grid3D> grid, const std::string &path,
                                                               WbWriter *const writer, SPtr<UbScheduler> s,
                                                               std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), comm(comm), writer(writer)
{
    init();
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensitySimulationObserver::init()
{
    gridRank     = grid->getRank();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);

        for (SPtr<Block3D> block : blockVector[level]) {
            UbTupleInt3 nx                           = grid->getBlockNX();
            SPtr<AverageValuesArray3D> averageValues = SPtr<AverageValuesArray3D>(
                new AverageValuesArray3D(val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 4, 0.0));
            block->getKernel()->getDataSet()->setAverageValues(averageValues);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensitySimulationObserver::update(real step)
{
    calculateAverageValues(int(step));

    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "TurbulenceIntensitySimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensitySimulationObserver::collectData(real step)
{
    int istep = int(step);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                addData(block);
            }
        }
    }

    std::string partName = writer->writeOctsWithNodeData(
        path + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep), nodes, cells, datanames, data);
    size_t found      = partName.find_last_of("//");
    std::string piece = partName.substr(found + 1);

    std::vector<std::string> cellDataNames;

    std::vector<std::string> pieces = comm->gather(piece);
    if (comm->getProcessID() == comm->getRoot()) {
        std::string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(
            path + "_" + UbSystem::toString(istep), pieces, datanames, cellDataNames);

        std::vector<std::string> filenames;
        filenames.push_back(pname);
        if (step == SimulationObserver::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(path + "_collection", filenames, istep, false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path + "_collection", filenames, istep, false);
        }
        UBLOG(logINFO, "TurbulenceIntensitySimulationObserver step: " << istep);
    }

    clearData();
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensitySimulationObserver::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensitySimulationObserver::addData(const SPtr<Block3D> block)
{
    UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
    //   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
    UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
    real dx                 = grid->getDeltaX(block);

    // Diese Daten werden geschrieben:
    datanames.resize(0);
    datanames.emplace_back("TI");

    data.resize(datanames.size());

    SPtr<ILBMKernel> kernel                 = block->getKernel();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageValues();
    // int ghostLayerWidth = kernel->getGhostLayerWidth();

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = int(distributions->getNX1());
    int maxX2 = int(distributions->getNX2());
    int maxX3 = int(distributions->getNX3());

    // nummern vergeben und node std::vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    // D3Q27BoundaryConditionPtr bcPtr;
    int nr = (int)nodes.size();

    for (int ix3 = minX3; ix3 < maxX3 - 1; ix3++) {
        for (int ix2 = minX2; ix2 < maxX2 - 1; ix2++) {
            for (int ix1 = minX1; ix1 < maxX1 - 1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    int index                  = 0;
                    nodeNumbers(ix1, ix2, ix3) = nr++;
                    nodes.push_back(makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1 * dx),
                                                float(val<2>(org) - val<2>(nodeOffset) + ix2 * dx),
                                                float(val<3>(org) - val<3>(nodeOffset) + ix3 * dx)));

                    // compute turbulence intensity
                    real temp =
                        (*av)(ix1, ix2, ix3, AvVxxyyzz) / ((*av)(ix1, ix2, ix3, AvVx) * (*av)(ix1, ix2, ix3, AvVx) +
                                                           (*av)(ix1, ix2, ix3, AvVy) * (*av)(ix1, ix2, ix3, AvVy) +
                                                           (*av)(ix1, ix2, ix3, AvVz) * (*av)(ix1, ix2, ix3, AvVz));

                    real ti = sqrt(temp);

                    if (UbMath::isNaN(ti))
                        UB_THROW(
                            UbException(UB_EXARGS, "TI is not a number (nan or -1.#IND), sqrt(temp), where temp = " +
                                                       UbSystem::toString(temp) +
                                                       ", AvVx = " + UbSystem::toString((*av)(ix1, ix2, ix3, AvVx)) +
                                                       " AvVy = " + UbSystem::toString((*av)(ix1, ix2, ix3, AvVy)) +
                                                       " AvVz = " + UbSystem::toString((*av)(ix1, ix2, ix3, AvVz))));

                    data[index++].push_back(ti);
                }
            }
        }
    }

    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

    // cell std::vector erstellen
    for (int ix3 = minX3; ix3 < maxX3 - 1; ix3++) {
        for (int ix2 = minX2; ix2 < maxX2 - 1; ix2++) {
            for (int ix1 = minX1; ix1 < maxX1 - 1; ix1++) {
                if ((SWB = nodeNumbers(ix1, ix2, ix3)) >= 0 && (SEB = nodeNumbers(ix1 + 1, ix2, ix3)) >= 0 &&
                    (NEB = nodeNumbers(ix1 + 1, ix2 + 1, ix3)) >= 0 && (NWB = nodeNumbers(ix1, ix2 + 1, ix3)) >= 0 &&
                    (SWT = nodeNumbers(ix1, ix2, ix3 + 1)) >= 0 && (SET = nodeNumbers(ix1 + 1, ix2, ix3 + 1)) >= 0 &&
                    (NET = nodeNumbers(ix1 + 1, ix2 + 1, ix3 + 1)) >= 0 &&
                    (NWT = nodeNumbers(ix1, ix2 + 1, ix3 + 1)) >= 0) {
                    cells.push_back(makeUbTuple((unsigned int)SWB, (unsigned int)SEB, (unsigned int)NEB,
                                                (unsigned int)NWB, (unsigned int)SWT, (unsigned int)SET,
                                                (unsigned int)NET, (unsigned int)NWT));
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensitySimulationObserver::calculateAverageValues(real timeStep)
{
    using namespace vf::lbm::dir;
    using namespace D3Q27System;

    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    real f[27];
    real vx, vy, vz;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                SPtr<ILBMKernel> kernel                 = block->getKernel();
                SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageValues();

                int minX1 = 0;
                int minX2 = 0;
                int minX3 = 0;

                int maxX1 = int(distributions->getNX1());
                int maxX2 = int(distributions->getNX2());
                int maxX3 = int(distributions->getNX3());

                for (int ix3 = minX3; ix3 < maxX3 - 1; ix3++) {
                    for (int ix2 = minX2; ix2 < maxX2 - 1; ix2++) {
                        for (int ix1 = minX1; ix1 < maxX1 - 1; ix1++) {
                            if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                                //////////////////////////////////////////////////////////////////////////
                                // read distribution
                                ////////////////////////////////////////////////////////////////////////////
                                distributions->getPreCollisionDistribution(f, ix1, ix2, ix3);
                                //////////////////////////////////////////////////////////////////////////
                                // compute velocity
                                //////////////////////////////////////////////////////////////////////////
                                vx = f[dP00] - f[dM00] + f[dPP0] - f[dMM0] + f[dPM0] - f[dMP0] + f[dP0P] - f[dM0M] + f[dP0M] - f[dM0P] +
                                     f[dPPP] - f[dMMP] + f[dPMP] - f[dMPP] + f[dPPM] - f[dMMM] + f[dPMM] - f[dMPM];

                                vy = f[d0P0] - f[d0M0] + f[dPP0] - f[dMM0] - f[dPM0] + f[dMP0] + f[d0PP] - f[d0MM] + f[d0PM] - f[d0MP] +
                                     f[dPPP] - f[dMMP] - f[dPMP] + f[dMPP] + f[dPPM] - f[dMMM] - f[dPMM] + f[dMPM];

                                vz = f[d00P] - f[d00M] + f[dP0P] - f[dM0M] - f[dP0M] + f[dM0P] + f[d0PP] - f[d0MM] - f[d0PM] + f[d0MP] +
                                     f[dPPP] + f[dMMP] + f[dPMP] + f[dMPP] - f[dPPM] - f[dMMM] - f[dPMM] - f[dMPM];
                                //////////////////////////////////////////////////////////////////////////
                                // compute average values
                                //////////////////////////////////////////////////////////////////////////
                                (*av)(ix1, ix2, ix3, AvVx) =
                                    ((*av)(ix1, ix2, ix3, AvVx) * timeStep + vx) / (timeStep + 1.0);
                                (*av)(ix1, ix2, ix3, AvVy) =
                                    ((*av)(ix1, ix2, ix3, AvVy) * timeStep + vy) / (timeStep + 1.0);
                                (*av)(ix1, ix2, ix3, AvVz) =
                                    ((*av)(ix1, ix2, ix3, AvVz) * timeStep + vz) / (timeStep + 1.0);

                                (*av)(ix1, ix2, ix3, AvVxxyyzz) =
                                    ((vx - (*av)(ix1, ix2, ix3, AvVx)) * (vx - (*av)(ix1, ix2, ix3, AvVx)) +
                                     (vy - (*av)(ix1, ix2, ix3, AvVy)) * (vy - (*av)(ix1, ix2, ix3, AvVy)) +
                                     (vz - (*av)(ix1, ix2, ix3, AvVz)) * (vz - (*av)(ix1, ix2, ix3, AvVz)) +
                                     (*av)(ix1, ix2, ix3, AvVxxyyzz) * timeStep) /
                                    (timeStep + 1.0);
                                //////////////////////////////////////////////////////////////////////////
                            }
                        }
                    }
                }
            }
        }
    }
}

//! \}
