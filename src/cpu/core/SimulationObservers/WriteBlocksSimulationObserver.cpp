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

#include "WriteBlocksSimulationObserver.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include <logger/Logger.h>

#include "Block3D.h"
#include <parallel/Communicator.h>
#include "D3Q27System.h"
#include "Grid3D.h"
#include "UbScheduler.h"

WriteBlocksSimulationObserver::WriteBlocksSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                               WbWriter *const writer, std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), writer(writer), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
WriteBlocksSimulationObserver::~WriteBlocksSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void WriteBlocksSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void WriteBlocksSimulationObserver::collectData(real step)
{
    if (comm->getProcessID() == comm->getRoot()) {
        int istep = int(step);
        std::vector<std::string> filenames;
        std::vector<UbTupleFloat3> nodes;
        std::vector<UbTupleInt8> cells;
        std::vector<std::string> celldatanames;

        celldatanames.emplace_back("isActive");
        celldatanames.emplace_back("rank");
        celldatanames.emplace_back("interface");
        celldatanames.emplace_back("ID");
        celldatanames.emplace_back("part");
        celldatanames.emplace_back("level");
        // celldatanames.push_back("connectorCF");
        // celldatanames.push_back("connectorFC");
#if defined VF_FETOL
        celldatanames.push_back("bundle");
#endif

        std::vector<std::vector<double>> celldata(celldatanames.size());

        int nr           = 0;
        int minInitLevel = this->grid->getCoarsestInitializedLevel();
        int maxInitLevel = this->grid->getFinestInitializedLevel();

        for (int level = minInitLevel; level <= maxInitLevel; level++) {
            std::vector<SPtr<Block3D>> blockVector;
            grid->getBlocks(level, blockVector);
            for (SPtr<Block3D> block : blockVector) {
                UbTupleDouble3 org         = grid->getBlockWorldCoordinates(block);
                UbTupleDouble3 blockLength = grid->getBlockLengths(block);

                nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org)), (float)(val<3>(org))));
                nodes.push_back(makeUbTuple((float)(val<1>(org) + val<1>(blockLength)), (float)(val<2>(org)),
                                            (float)(val<3>(org))));
                nodes.push_back(makeUbTuple((float)(val<1>(org) + val<1>(blockLength)),
                                            (float)(val<2>(org) + val<2>(blockLength)), (float)(val<3>(org))));
                nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org) + val<2>(blockLength)),
                                            (float)(val<3>(org))));
                nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org)),
                                            (float)(val<3>(org) + val<3>(blockLength))));
                nodes.push_back(makeUbTuple((float)(val<1>(org) + val<1>(blockLength)), (float)(val<2>(org)),
                                            (float)(val<3>(org) + val<3>(blockLength))));
                nodes.push_back(makeUbTuple((float)(val<1>(org) + val<1>(blockLength)),
                                            (float)(val<2>(org) + val<2>(blockLength)),
                                            (float)(val<3>(org) + val<3>(blockLength))));
                nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org) + val<2>(blockLength)),
                                            (float)(val<3>(org) + val<3>(blockLength))));
                cells.push_back(makeUbTuple(nr, nr + 1, nr + 2, nr + 3, nr + 4, nr + 5, nr + 6, nr + 7));
                nr += 8;

                // data
                celldata[0].push_back((real)block->isActive());
                celldata[1].push_back((real)block->getRank());
                celldata[2].push_back((real)block->hasInterpolationFlag());
                celldata[3].push_back((real)block->getGlobalID());
                celldata[4].push_back((real)block->getPart());
                celldata[5].push_back((real)block->getLevel());

                // bool flag = false;
                // std::vector<SPtr<Block3DConnector>> connectors;

                // block->pushBackLocalInterpolationConnectorsCF(connectors);
                // for (std::size_t i = 0; i<connectors.size(); i++)
                //   if (connectors[i])
                //   {
                //      if (connectors[i]->getSendDir() == D3Q27System::d0MM)
                //      {

                //         flag = true;
                //      }
                //   }

                // if (flag)
                //{
                //   celldata[6].push_back(1);
                //   UBLOG(logINFO, "CF: "+block->toString());
                //}
                // else
                //{
                //   celldata[6].push_back(0);
                //}

                // flag = false;
                // connectors.resize(0);
                // block->pushBackLocalInterpolationConnectorsFC(connectors);
                // for (std::size_t i = 0; i<connectors.size(); i++)
                //   if (connectors[i])
                //   {
                //      if (connectors[i]->getSendDir() == D3Q27System::d0MM)
                //      {

                //         flag = true;
                //      }
                //   }

                // if (flag)
                //{
                //   celldata[7].push_back(1);
                //   UBLOG(logINFO, "FC: "+block->toString());
                //}
                // else
                //{
                //   celldata[7].push_back(0);
                //}

#ifdef VF_FETOL
                celldata[6].push_back((real)block->getBundle());
#endif
            }
        }

        filenames.push_back(writer->writeOctsWithCellData(
            path + "/blocks/blocks_" + UbSystem::toString(grid->getRank()) + "_" + UbSystem::toString(istep), nodes,
            cells, celldatanames, celldata));

        if (istep == SimulationObserver::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(path + "/blocks/blocks_collection", filenames, istep,
                                                                false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path + "/blocks/blocks_collection", filenames,
                                                                     istep, false);
        }

        VF_LOG_INFO("WriteBlocksSimulationObserver step: {}", istep);
    }
}

//! \}
