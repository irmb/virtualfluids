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
#ifndef QVTKWRITER_HPP
#define QVTKWRITER_HPP

#include <array>
#include <vector>

#include <basics/StringUtilities/StringUtil.h>
#include <basics/utilities/UbSystem.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <lbm/constants/D3Q27.h>
#include <logger/Logger.h>

#include "gpu/GridGenerator/grid/NodeValues.h"
#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Utilities/FindNeighbors.h"

namespace QDebugVtkWriter
{

using namespace vf::lbm::dir;

namespace
{
inline void modifyLineLengthsForQs(const std::array<double, 3> &coords, std::array<double, 3> &neighborCoords, real q)
{
    if (q == 1.0 || q <= 0.0)
        return;

    const auto dx = neighborCoords[0] - coords[0];
    const auto dy = neighborCoords[1] - coords[1];
    const auto dz = neighborCoords[2] - coords[2];

    neighborCoords[0] = coords[0] + q * dx;
    neighborCoords[1] = coords[1] + q * dy;
    neighborCoords[2] = coords[2] + q * dz;
}

inline void writeQLines(LBMSimulationParameter *parH, QforBoundaryConditions &boundaryQ, const std::string &filepath,
                        WbWriter *writer)
{
    VF_LOG_INFO("Write qs in for boundary condition to {}.", filepath);

    const auto numberOfNodes = boundaryQ.numberOfBCnodes;
    std::vector<UbTupleFloat3> nodes;
    nodes.reserve(numberOfNodes * 8 * 2);
    std::vector<UbTupleInt2> lines;
    lines.reserve(numberOfNodes * 8);

    std::vector<std::string> dataNames = { "nodeIndex", "q" };
    std::vector<std::vector<float>> lineData(2);

    for (size_t i = 0; i < numberOfNodes; i++) {
        const auto nodeIndex = boundaryQ.k[i];
        const std::array<double, 3> coords = { parH->coordinateX[nodeIndex], parH->coordinateY[nodeIndex],
                                               parH->coordinateZ[nodeIndex] };

        for (size_t direction = 1; direction < ENDDIR; direction++) {

            const auto q = boundaryQ.q27[direction][i];
            if (q <= (real)0.0) {
                continue;
            }

            const auto positionNeighbor = getNeighborIndex(parH, (uint)nodeIndex, (int)direction);

            std::array<double, 3> neighborCoords = { parH->coordinateX[positionNeighbor],
                                                     parH->coordinateY[positionNeighbor],
                                                     parH->coordinateZ[positionNeighbor] };

            modifyLineLengthsForQs(coords, neighborCoords, q);

            nodes.emplace_back(float(coords[0]), float(coords[1]), coords[2]);
            nodes.emplace_back(float(neighborCoords[0]), float(neighborCoords[1]), float(neighborCoords[2]));

            lines.emplace_back((int)nodes.size() - 2, (int)nodes.size() - 1);
            lineData[0].push_back(nodeIndex);
            lineData[1].push_back(q);
        }
    }

    writer->writeLinesWithLineData(filepath, nodes, lines, dataNames, lineData);
}
} // namespace

inline void writeQLinesDebug(Parameter *para, QforBoundaryConditions &boundaryQ, uint level, const std::string& fileName)
{
    const auto filePath = para->getFName() + "_" + fileName + ".vtk";
    auto writer = WbWriterVtkXmlBinary::getInstance();
    writeQLines(para->getParH(level).get(), boundaryQ, filePath, writer);
}

} // namespace QDebugVtkWriter

#endif

//! \}
