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

#include "WriteMacroscopicQuantitiesSimulationObserver.h"
#include "BCSet.h"
#include "LBMKernel.h"
#include <string>
#include <vector>

#include "BCArray3D.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "DataSet3D.h"
#include "Grid3D.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

WriteMacroscopicQuantitiesSimulationObserver::WriteMacroscopicQuantitiesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
WriteMacroscopicQuantitiesSimulationObserver::WriteMacroscopicQuantitiesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                                             const std::string &path,
                                                                             WbWriter *const writer,
                                                                             SPtr<LBMUnitConverter> conv,
                                                                             std::shared_ptr<vf::parallel::Communicator> comm)
        : SimulationObserver(grid, s), path(path), writer(writer), conv(conv), comm(comm)
{
    gridRank = comm->getProcessID();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++)
    {
        grid->getBlocks(level, gridRank, true, blockVector[level]);
    }
}

//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesSimulationObserver::init()
{}

//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "WriteMacroscopicQuantitiesSimulationObserver::update:" << step);
}

//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesSimulationObserver::collectData(real step)
{
    int istep = static_cast<int>(step);

    for (int level = minInitLevel; level <= maxInitLevel; level++)
    {
        for (SPtr<Block3D> block : blockVector[level])
        {
            if (block)
            {
                addDataMQ(block);
            }
        }
    }

    std::string pfilePath, partPath, subfolder, cfilePath;

    subfolder = "mq" + UbSystem::toString(istep);
    pfilePath = path + "/mq/" + subfolder;
    cfilePath = path + "/mq/mq_collection";
    partPath = pfilePath + "/mq" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

    std::string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
    size_t found = partName.find_last_of("/");
    std::string piece = partName.substr(found + 1);
    piece = subfolder + "/" + piece;

    std::vector<std::string> cellDataNames;
    std::vector<std::string> pieces = comm->gather(piece);
    if (comm->getProcessID() == comm->getRoot()) {
        std::string pname =
                WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
        found = pname.find_last_of("/");
        piece = pname.substr(found + 1);

        std::vector<std::string> filenames;
        filenames.push_back(piece);
        if (step == SimulationObserver::scheduler->getMinBegin())
        {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else
        {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        UBLOG(logINFO, "WriteMacroscopicQuantitiesSimulationObserver step: " << istep);
    }

    clearData();
}

//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesSimulationObserver::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}

//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesSimulationObserver::addDataMQ(SPtr<Block3D> block)
{
    real level   = (real)block->getLevel();

    // This data is written:
    datanames.resize(0);
    datanames.push_back("Rho");
    datanames.push_back("Vx");
    datanames.push_back("Vy");
    datanames.push_back("Vz");
    datanames.push_back("Level");

    data.resize(datanames.size());

    SPtr<ILBMKernel> kernel                 = block->getKernel();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    real f[D3Q27System::ENDF + 1];
    real vx1, vx2, vx3, rho;

    // Node numbering always starts at 0!
    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

    if (block->getKernel()->getCompressible()) {
        calcMacros = &D3Q27System::calcCompMacroscopicValues;
    } else {
        calcMacros = &D3Q27System::calcIncompMacroscopicValues;
    }

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = (int)(distributions->getNX1());
    int maxX2 = (int)(distributions->getNX2());
    int maxX3 = (int)(distributions->getNX3());

    // Assign numbers and create node vector + collect data
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    maxX1 -= 2;
    maxX2 -= 2;
    maxX3 -= 2;

    int nr = (int)nodes.size();

    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    int index                  = 0;
                    nodeNumbers(ix1, ix2, ix3) = nr++;
                    Vector3D worldCoordinates  = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                    nodes.push_back(UbTupleFloat3(real(worldCoordinates[0]), real(worldCoordinates[1]),
                                                  real(worldCoordinates[2])));

                    distributions->getPreCollisionDistribution(f, ix1, ix2, ix3);
                    calcMacros(f, rho, vx1, vx2, vx3);

                    if (UbMath::isNaN(rho) || UbMath::isInfinity(rho))
                         UB_THROW( UbException(UB_EXARGS,"rho is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+",node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                    if (UbMath::isNaN(vx1) || UbMath::isInfinity(vx1))
                         UB_THROW( UbException(UB_EXARGS,"vx1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                    if (UbMath::isNaN(vx2) || UbMath::isInfinity(vx2))
                         UB_THROW( UbException(UB_EXARGS,"vx2 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                    if (UbMath::isNaN(vx3) || UbMath::isInfinity(vx3))
                         UB_THROW( UbException(UB_EXARGS,"vx3 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));

                    data[index++].push_back(rho);
                    data[index++].push_back(vx1);
                    data[index++].push_back(vx2);
                    data[index++].push_back(vx3);

                    data[index++].push_back(level);
                }
            }
        }
    }
    maxX1 -= 1;
    maxX2 -= 1;
    maxX3 -= 1;
    // create cell vector 
    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
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

//! \}
