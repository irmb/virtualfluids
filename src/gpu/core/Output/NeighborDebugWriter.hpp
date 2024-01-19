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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Anna Wellmann
//=======================================================================================
#ifndef NEIGHBORDEBUG_HPP
#define NEIGHBORDEBUG_HPP

#include <logger/Logger.h>

#include <basics/utilities/UbSystem.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <GridGenerator/grid/NodeValues.h>

#include <lbm/constants/D3Q27.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include "StringUtilities/StringUtil.h"
#include "Utilities/FindNeighbors.h"

namespace NeighborDebugWriter
{

inline void writeNeighborLinkLines(LBMSimulationParameter *parH, int direction, const std::string &name, WbWriter *writer)
{
    VF_LOG_INFO("Write node links in direction {}.", direction);

    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2> cells;

    for (size_t position = 0; position < parH->numberOfNodes; position++) {
        if (parH->typeOfGridNode[position] != GEO_FLUID) continue;

        const double x1 = parH->coordinateX[position];
        const double x2 = parH->coordinateY[position];
        const double x3 = parH->coordinateZ[position];

        const uint positionNeighbor = getNeighborIndex(parH, (uint)position, direction);

        const double x1Neighbor = parH->coordinateX[positionNeighbor];
        const double x2Neighbor = parH->coordinateY[positionNeighbor];
        const double x3Neighbor = parH->coordinateZ[positionNeighbor];

        nodes.emplace_back(float(x1), float(x2), float(x3));
        nodes.emplace_back(float(x1Neighbor), float(x2Neighbor), float(x3Neighbor));

        cells.emplace_back((int)nodes.size() - 2, (int)nodes.size() - 1);
    }
    writer->writeLines(name, nodes, cells);
}

inline void writeNeighborLinkLinesDebug(Parameter *para)
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (size_t direction = vf::lbm::dir::STARTDIR; direction <= vf::lbm::dir::ENDDIR; direction++) {
            const std::string fileName = para->getFName() + "_" + StringUtil::toString<int>(level) + "_Link_" +
                                         std::to_string(direction) + "_Debug.vtk";
            writeNeighborLinkLines(para->getParH(level).get(), (int)direction, fileName,
                                   WbWriterVtkXmlBinary::getInstance());
        }
    }
}

} // namespace NeighborDebugWriter

#endif

//! \}
