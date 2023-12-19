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

#include "WriteBoundaryConditionsSimulationObserver.h"
#include "BCSet.h"
#include "LBMKernel.h"
#include <string>
#include <vector>

#include <logger/Logger.h>

#include "BCArray3D.h"
#include "Block3D.h"
#include "CbArray3D.h"
#include <parallel/Communicator.h>
#include "Grid3D.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"
#include "WbWriter.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

using namespace std;

WriteBoundaryConditionsSimulationObserver::WriteBoundaryConditionsSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
WriteBoundaryConditionsSimulationObserver::WriteBoundaryConditionsSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                                       const std::string &path, WbWriter *const writer,
                                                                       std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), writer(writer), comm(comm)
{
    gridRank     = comm->getProcessID();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);
    }
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "WriteBoundaryConditionsSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsSimulationObserver::collectData(real step)
{
    int istep = static_cast<int>(step);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                addDataGeo(block);
            }
        }
    }

    string pfilePath, partPath, subfolder, cfilePath;

    subfolder = "bc" + UbSystem::toString(istep);
    pfilePath = path + "/bc/" + subfolder;
    cfilePath = path + "/bc/bc_collection";
    partPath  = pfilePath + "/bc" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

    string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
    size_t found    = partName.find_last_of("/");
    string piece    = partName.substr(found + 1);
    piece           = subfolder + "/" + piece;

    vector<string> cellDataNames;
    vector<string> pieces = comm->gather(piece);
    if (comm->getProcessID() == comm->getRoot()) {
        string pname =
            WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
        found = pname.find_last_of("/");
        piece = pname.substr(found + 1);

        vector<string> filenames;
        filenames.push_back(piece);
        if (step == SimulationObserver::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        VF_LOG_INFO("WriteBoundaryConditionsSimulationObserver step: {}", istep);
    }

    clearData();
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsSimulationObserver::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsSimulationObserver::addDataGeo(SPtr<Block3D> block)
{
    using namespace vf::basics::constant;

    UbTupleDouble3 org        = grid->getBlockWorldCoordinates(block);
    UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
    real dx                 = grid->getDeltaX(block);

    real level = (real)block->getLevel();

    // Diese Daten werden geschrieben:
    datanames.resize(0);
    datanames.emplace_back("Boundary Conditions");
    datanames.emplace_back("Geometry");
    datanames.emplace_back("Level");

    data.resize(datanames.size());

    SPtr<ILBMKernel> kernel = block->getKernel();
    SPtr<BCArray3D> bcArray = kernel->getBCSet()->getBCArray();


    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = (int)bcArray->getNX1();
    int maxX2 = (int)bcArray->getNX2();
    int maxX3 = (int)bcArray->getNX3();

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    // D3Q27BoundaryConditionPtr bcPtr;
    int nr = (int)nodes.size();

    maxX1 -= 1;
    maxX2 -= 1;
    maxX3 -= 1;

    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3)) {

                    // int index = 0;
                    nodeNumbers(ix1, ix2, ix3) = nr++;
                    nodes.push_back(makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1 * dx),
                                                float(val<2>(org) - val<2>(nodeOffset) + ix2 * dx),
                                                float(val<3>(org) - val<3>(nodeOffset) + ix3 * dx)));

                    auto bc = bcArray->getBC(ix1, ix2, ix3);
                    if (!bcArray->hasBC(ix1, ix2, ix3)) {
                        data[0].push_back(c0o1);
                    } else if (bc && bc->hasNoSlipBoundary())
                        data[0].push_back(c1o1);
                    else if (bc && bc->hasVelocityBoundary())
                        data[0].push_back(c2o1);
                    else if (bc && bc->hasDensityBoundary())
                        data[0].push_back(c3o1);
                    else if (bc && bc->hasSlipBoundary())
                        data[0].push_back(c4o1);
                    // else
                    //   data[0].push_back(5.0);

                    if (bcArray->isSolid(ix1, ix2, ix3)) {
                        data[1].push_back(c1o1);
                    } else {
                        data[1].push_back(c0o1);
                    }

                    data[2].push_back(level);
                }
            }
        }
    }

    maxX1 -= 1;
    maxX2 -= 1;
    maxX3 -= 1;

    // knotennummerierung faengt immer bei 0 an!
    int SWB = 0, SEB = 0, NEB = 0, NWB = 0, SWT = 0, SET = 0, NET = 0, NWT = 0;

    // cell vector erstellen
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
