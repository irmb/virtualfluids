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
//! \author Hussein Alihussein
//! \note Grid-winding diagnostics IO helpers.
//=======================================================================================
#pragma once

#include <basics/PointerDefinitions.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <logger/Logger.h>

#include <parallel/Communicator.h>

#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace vf::grid_winding
{

//==============================================================================
// Missing-link diagnostics (common CPU/GPU)
//==============================================================================

inline MissingLinkCollection collectMissingLinks(const SubgridDistanceStats&             stats,
                                                 const SPtr<vf::parallel::Communicator>& comm)
{
    MissingLinkCollection collection;

    const int localCount = static_cast<int>(stats.missingLinks.size());
    std::vector<int> sendCounts{ localCount };
    std::vector<int> recvCounts;

    if (comm)
        comm->allGather(sendCounts, recvCounts);

    if (recvCounts.empty())
        recvCounts = sendCounts;

    std::vector<double> coordsSend;
    coordsSend.reserve(static_cast<std::size_t>(localCount) * 6);
    std::vector<int>   dirsSend;
    dirsSend.reserve(static_cast<std::size_t>(localCount));
    std::vector<float> windSend;
    windSend.reserve(static_cast<std::size_t>(localCount));
    std::vector<float> endpointSend;
    endpointSend.reserve(static_cast<std::size_t>(localCount));

    for (const auto& missing : stats.missingLinks) {
        coordsSend.push_back(missing.origin[0]);
        coordsSend.push_back(missing.origin[1]);
        coordsSend.push_back(missing.origin[2]);
        coordsSend.push_back(missing.neighbour[0]);
        coordsSend.push_back(missing.neighbour[1]);
        coordsSend.push_back(missing.neighbour[2]);
        dirsSend.push_back(missing.direction);
        windSend.push_back(missing.windingNumber);
        endpointSend.push_back(missing.endpointKind);
    }

    std::vector<double> coordsRecv;
    std::vector<int>    dirsRecv;
    std::vector<float>  windRecv;
    std::vector<float>  endpointRecv;
    if (comm) {
        comm->allGather(coordsSend, coordsRecv);
        comm->allGather(dirsSend, dirsRecv);
        comm->allGather(windSend, windRecv);
        comm->allGather(endpointSend, endpointRecv);
    }

    if (coordsRecv.empty())
        coordsRecv = coordsSend;
    if (dirsRecv.empty())
        dirsRecv = dirsSend;
    if (windRecv.empty())
        windRecv = windSend;
    if (endpointRecv.empty())
        endpointRecv = endpointSend;

    const std::size_t totalMissing = std::accumulate(recvCounts.begin(), recvCounts.end(), std::size_t{ 0 });

    if (totalMissing == 0)
        return collection;

    if (comm && !comm->isRoot())
        return collection;

    const std::size_t gatheredCount = dirsRecv.size();
    if (gatheredCount == 0)
        return collection;

    if (coordsRecv.size() < gatheredCount * 6) {
        UBLOG(logWARNING, "collectMissingLinks: received coordinate buffer too small (" << coordsRecv.size()
                                                                                        << " values for "
                                                                                        << gatheredCount
                                                                                        << " links)");
        return collection;
    }

    collection.nodes.reserve(gatheredCount * 2);
    collection.lines.reserve(gatheredCount);
    collection.directions.reserve(gatheredCount);
    collection.windingNumbers.reserve(gatheredCount);
    collection.endpointKind.reserve(gatheredCount);

    for (std::size_t idx = 0; idx < gatheredCount; ++idx) {
        const std::size_t base = idx * 6;
        const float       ox   = static_cast<float>(coordsRecv[base + 0]);
        const float       oy   = static_cast<float>(coordsRecv[base + 1]);
        const float       oz   = static_cast<float>(coordsRecv[base + 2]);
        const float       nx   = static_cast<float>(coordsRecv[base + 3]);
        const float       ny   = static_cast<float>(coordsRecv[base + 4]);
        const float       nz   = static_cast<float>(coordsRecv[base + 5]);

        const int nodeIndex = static_cast<int>(collection.nodes.size());
        collection.nodes.emplace_back(ox, oy, oz);
        collection.nodes.emplace_back(nx, ny, nz);
        collection.lines.emplace_back(nodeIndex, nodeIndex + 1);
        collection.directions.push_back(static_cast<float>(dirsRecv[idx]));
        const float wind = idx < windRecv.size() ? windRecv[idx]
                                                                  : std::numeric_limits<float>::quiet_NaN();
        collection.windingNumbers.push_back(wind);
        const float endpoint = idx < endpointRecv.size() ? endpointRecv[idx] : 0.0f;
        collection.endpointKind.push_back(endpoint);
    }

    collection.total = collection.lines.size();
    return collection;
}

inline void populateMissingLinkEndpointKind(MissingLinkCollection &collection, const SPtr<GbTriFaceMesh3D> &surface)
{
    if (collection.total == 0 || !surface)
        return;

    collection.endpointKind.assign(collection.total, 0.0f);

    for (std::size_t idx = 0; idx < collection.total; ++idx) {
        const auto &originNode    = collection.nodes[idx * 2 + 0];
        const auto &neighbourNode = collection.nodes[idx * 2 + 1];

        Vec3 origin(static_cast<double>(val<1>(originNode)),
                    static_cast<double>(val<2>(originNode)),
                    static_cast<double>(val<3>(originNode)));
        Vec3 neighbour(static_cast<double>(val<1>(neighbourNode)),
                       static_cast<double>(val<2>(neighbourNode)),
                       static_cast<double>(val<3>(neighbourNode)));

        const double maxDistance = (neighbour - origin).Length();
        const double eps = std::max(1.0e-4, 1.0e-3 * maxDistance);

        const int kind = vf::grid_winding::endpointSurfaceKind(*surface, origin, neighbour, eps);
        collection.endpointKind[idx] = static_cast<float>(kind);
    }
}

inline void writeMissingLinks(const MissingLinkCollection& collection, const std::string& basePath)
{
    if (collection.empty())
        return;

    auto nodes = collection.nodes;
    auto lines = collection.lines;

    std::vector<std::vector<float>> lineData;
    lineData.reserve(3);
    lineData.emplace_back(collection.directions);

    std::vector<std::string> lineDataNames;
    lineDataNames.reserve(3);
    lineDataNames.emplace_back("lattice_dir");

    if (collection.endpointKind.size() == collection.total) {
        lineData.emplace_back(collection.endpointKind);
        lineDataNames.emplace_back("endpoint_kind");
    }

    if (collection.windingNumbers.size() == collection.total) {
        lineData.emplace_back(collection.windingNumbers);
        lineDataNames.emplace_back("winding");
    }

    const std::string filename = WbWriterVtkXmlBinary::getInstance()->writeLinesWithLineData(
        basePath, nodes, lines, lineDataNames, lineData);
    UBLOG(logINFO, "Wrote " << collection.total << " missing boundary links to " << filename);
}

//==============================================================================
// Q-line diagnostics (common CPU/GPU)
//==============================================================================

//! \brief Build a QLineCollection from gathered coordinate / direction / q buffers.
//!
//! This helper is intentionally platform agnostic; CPU and GPU code are responsible
//! for filling the input buffers and performing any MPI all-gather before calling it.
inline QLineCollection buildQLineCollection(const std::vector<double>& coordsRecv,
                                            const std::vector<int>&    dirRecv,
                                            const std::vector<float>&  qRecv)
{
    QLineCollection collection;

    if (coordsRecv.empty() || dirRecv.empty() || qRecv.empty())
        return collection;

    const std::size_t lineCount = qRecv.size();
    if (coordsRecv.size() < lineCount * 6)
        return collection;

    collection.nodes.reserve(lineCount * 2);
    collection.lines.reserve(lineCount);
    collection.directions.reserve(lineCount);
    collection.qValues.reserve(lineCount);
    collection.surfacePoints.reserve(lineCount);

    for (std::size_t idx = 0; idx < lineCount; ++idx) {
        const std::size_t base = idx * 6;

        const float ox = static_cast<float>(coordsRecv[base + 0]);
        const float oy = static_cast<float>(coordsRecv[base + 1]);
        const float oz = static_cast<float>(coordsRecv[base + 2]);
        const float hx = static_cast<float>(coordsRecv[base + 3]);
        const float hy = static_cast<float>(coordsRecv[base + 4]);
        const float hz = static_cast<float>(coordsRecv[base + 5]);

        const int startIndex = static_cast<int>(collection.nodes.size());
        collection.nodes.emplace_back(ox, oy, oz);
        collection.nodes.emplace_back(hx, hy, hz);
        collection.lines.emplace_back(startIndex, startIndex + 1);

        const int dir = idx < dirRecv.size() ? dirRecv[idx] : 0;
        collection.directions.push_back(static_cast<float>(dir));
        collection.qValues.push_back(qRecv[idx]);
        collection.surfacePoints.emplace_back(hx, hy, hz);
    }

    collection.total = collection.lines.size();
    return collection;
}

inline void writeQLines(const QLineCollection& collection, const std::string& basePath)
{
    if (collection.empty())
        return;

    const std::size_t lineCount = collection.lines.size();
    VF_LOG_INFO("[q-export] lines(q>0)={}", lineCount);

    auto nodes = collection.nodes;
    auto lines = collection.lines;

    std::vector<std::vector<float>> lineData(2);
    lineData[0] = collection.directions;
    lineData[1] = collection.qValues;
    std::vector<std::string> lineDataNames{ "direction", "q" };

    const std::string filename = WbWriterVtkXmlBinary::getInstance()->writeLinesWithLineData(
        basePath, nodes, lines, lineDataNames, lineData);
    UBLOG(logINFO, "Wrote " << collection.total << " boundary q lines to " << filename);

    if (!collection.surfacePoints.empty()) {
        auto pointNodes = collection.surfacePoints;
        std::vector<std::vector<double>> pointData(1);
        pointData[0].assign(collection.qValues.begin(), collection.qValues.end());
        std::vector<std::string> dataNames{ "q" };
        const std::string pointFilename = WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(
            basePath + "_surface_points", pointNodes, dataNames, pointData);
        UBLOG(logINFO, "Wrote " << pointNodes.size() << " surface points to " << pointFilename);
    }

}

} // namespace vf::grid_winding
