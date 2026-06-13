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
//! \brief Helper functions for GPU boundary-link processing and Q computation on triangle surfaces.
//! \note Generated with the assistance of OpenAI Codex (GPT-5). Reviewed and adapted by author.
//=======================================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <basics/Timer/Timer.h>
#include <logger/Logger.h>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <geometry3d/GbTriFaceMesh3D.h>

#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>

#include "grid/Grid.h"
#include "grid/GridImp.h"
#include "grid/NodeValues.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace vf::gpu::grid_winding
{
namespace detail
{
inline bool isSolidNeighbour(char neighbourType)
{
    switch (neighbourType) {
    case vf::gpu::STOPPER_SOLID:
    case vf::gpu::INVALID_SOLID:
        return true;
    default:
        return false;
    }
}
} // namespace detail

struct BoundaryProcessingConfig
{
    bool rebuildBoundaryQIndices{ true };
    bool ensureQStorage{ true };
    bool fillMissingQs{ true };
    real missingQValue{ static_cast<real>(vf::grid_winding::defaultMissingQ()) };
    bool parallelize{ true };
    std::size_t parallelThreshold{ 4096 };
};

inline void processSolidBoundaryLinks(GridImp& grid,
                                      const vf::grid_winding::TriangleSurfaces& surfaces,
                                      vf::grid_winding::SubgridDistanceStats* stats = nullptr)
{
    using vf::grid_winding::SubgridDistanceStats;
    using vf::grid_winding::Vec3;

    if (surfaces.empty())
        return;

    const double dx       = static_cast<double>(grid.getDelta());
    const int    dirStart = grid.getStartDirection();
    const int    dirEnd   = grid.getEndDirection();

    const uint nodeCount = grid.getSize();
    for (uint index = 0; index < nodeCount; ++index) {
        if (grid.getFieldEntry(index) != vf::gpu::BC_SOLID)
            continue;

        real ox = 0.0;
        real oy = 0.0;
        real oz = 0.0;
        grid.transIndexToCoords(index, ox, oy, oz);

        Vec3 origin{ static_cast<double>(ox), static_cast<double>(oy), static_cast<double>(oz) };

        for (int dir = dirStart; dir <= dirEnd; ++dir) {
            const auto& direction = grid.distribution.directions[dir];
            const int   cx        = direction[0];
            const int   cy        = direction[1];
            const int   cz        = direction[2];
            if (cx == 0 && cy == 0 && cz == 0)
                continue;

            const double neighbourX = origin.X1() + dx * static_cast<double>(cx);
            const double neighbourY = origin.X2() + dx * static_cast<double>(cy);
            const double neighbourZ = origin.X3() + dx * static_cast<double>(cz);

            const uint neighbourIndex = grid.transCoordToIndex(static_cast<real>(neighbourX),
                                                               static_cast<real>(neighbourY),
                                                               static_cast<real>(neighbourZ));
            if (neighbourIndex == INVALID_INDEX)
                continue;

            const char neighbourType = grid.getFieldEntry(neighbourIndex);
            if (!detail::isSolidNeighbour(neighbourType))
                continue;

            if (stats)
                ++stats->totalLinks;

            auto promoteSolid = [&]() -> bool {
                if (grid.getFieldEntry(index) == vf::gpu::FLUID) {
                    grid.setFieldEntry(index, vf::gpu::BC_SOLID);
                    return true;
                }
                return false;
            };

            auto existingQ = [&]() -> double {
                if (!grid.hasQIndex(index))
                    return -1.0;
                return static_cast<double>(grid.getQValue(index, static_cast<uint>(dir)));
            };

            auto setQ = [&](double q, int patchIndex) {
                if (!grid.hasQIndex(index))
                    return;
                grid.setQValue(index, dir, static_cast<real>(q));
                if (patchIndex >= 0)
                    grid.setQPatch(index, static_cast<uint>(patchIndex));
            };

            auto markSurface = [](double, int) {};

            auto getWindingValue = []() -> std::pair<bool, double> {
                return { false, 0.0 };
            };

            auto handleMissing = [&](int missingDir, const Vec3&, const Vec3&,
                                     SubgridDistanceStats::MissingLink&) -> bool {
                if (!grid.hasQIndex(index))
                    return false;

                const auto& directionVec = grid.distribution.directions[missingDir];
                const int   cxMiss       = directionVec[0];
                const int   cyMiss       = directionVec[1];
                const int   czMiss       = directionVec[2];
                if (cxMiss == 0 && cyMiss == 0 && czMiss == 0)
                    return false;

                real nx = static_cast<real>(origin.X1()) + static_cast<real>(cxMiss) * grid.getDelta();
                real ny = static_cast<real>(origin.X2()) + static_cast<real>(cyMiss) * grid.getDelta();
                real nz = static_cast<real>(origin.X3()) + static_cast<real>(czMiss) * grid.getDelta();
                const uint neighbourIdx = grid.transCoordToIndex(nx, ny, nz);
                if (neighbourIdx == INVALID_INDEX)
                    return false;

                const char neighbourTypeMiss = grid.getFieldEntry(neighbourIdx);
                if (!detail::isSolidNeighbour(neighbourTypeMiss))
                    return false;

                const double fallbackQ = vf::grid_winding::defaultMissingQ();
                grid.setQValue(index, missingDir, static_cast<real>(fallbackQ));
                return true;
            };

            vf::grid_winding::processBoundaryLink(surfaces,
                                                  origin,
                                                  dx,
                                                  cx,
                                                  cy,
                                                  cz,
                                                  dir,
                                                  promoteSolid,
                                                  existingQ,
                                                  setQ,
                                                  markSurface,
                                                  getWindingValue,
                                                  handleMissing,
                                                  stats);
        }
    }
}

inline void mergeSubgridStats(vf::grid_winding::SubgridDistanceStats &dst,
                              const vf::grid_winding::SubgridDistanceStats &src)
{
    dst.totalLinks += src.totalLinks;
    dst.hitLinks += src.hitLinks;
    dst.promotedFluidNodes += src.promotedFluidNodes;
    dst.promotedZeroQNodes += src.promotedZeroQNodes;
    dst.skippedNonLocalLinks += src.skippedNonLocalLinks;
    dst.missingLinks.insert(dst.missingLinks.end(), src.missingLinks.begin(), src.missingLinks.end());
}

inline void processSolidBoundaryLinksMultithreaded(GridImp &grid,
                                                   const vf::grid_winding::TriangleSurfaces& surfaces,
                                                   vf::grid_winding::SubgridDistanceStats *stats,
                                                   bool parallelize,
                                                   std::size_t parallelThreshold,
                                                   const std::vector<GbTriFaceMesh3D *> *meshes = nullptr)
{
    if (surfaces.empty())
        return;

    const uint nodeCount = grid.getSize();

    if (!parallelize || nodeCount < parallelThreshold) {
        vf::gpu::grid_winding::processSolidBoundaryLinks(grid, surfaces, stats);
        return;
    }

    vf::basics::Timer boundaryPrepTimer;
    boundaryPrepTimer.start();

    const uint boundaryCount = grid.getNumberOfSolidBoundaryNodes();
    if (boundaryCount == 0)
        return;

    std::vector<kd_tree::Tree<double> *> kdTrees;
    if (meshes && !meshes->empty())
        kdTrees = vf::grid_winding::ensureKdTrees(*meshes);

    if (meshes && !meshes->empty() && kdTrees.empty()) {
        VF_LOG_WARNING("No kd-tree available for grid-winding; boundary links will be treated as misses.");
    } else {
        VF_LOG_TRACE("Boundary processing using {} kd-trees", kdTrees.size());
    }

    VF_LOG_TRACE("Boundary indexing pass started: nodeCount={}, solidBoundaryNodes={}",
                 nodeCount, boundaryCount);
    std::vector<uint> boundaryIndices(boundaryCount, std::numeric_limits<uint>::max());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int idx = 0; idx < static_cast<int>(nodeCount); ++idx) {
        const uint index = static_cast<uint>(idx);
        if (!grid.hasQIndex(index))
            continue;
        const uint qIndex = grid.getQIndex(index);
        if (qIndex >= boundaryCount)
            continue;
        boundaryIndices[qIndex] = index;
    }
    VF_LOG_TRACE("Boundary indexing pass populated {} entries (expected {})",
                 boundaryIndices.size(), boundaryCount);

    boundaryIndices.erase(std::remove(boundaryIndices.begin(), boundaryIndices.end(), std::numeric_limits<uint>::max()),
                          boundaryIndices.end());

    const double prepSeconds = boundaryPrepTimer.getCurrentRuntimeInSeconds();
    VF_LOG_TRACE("Boundary indexing pass populated {} valid entries in {} sec", boundaryIndices.size(), prepSeconds);

    if (boundaryIndices.empty())
        return;

#if defined(_OPENMP)
    omp_set_dynamic(0);
    vf::basics::Timer boundaryTimer;
    boundaryTimer.start();
    const double dx          = static_cast<double>(grid.getDelta());
    const int    dirStart    = grid.getStartDirection();
    const int    dirEnd      = grid.getEndDirection();
    const bool   collectStats = (stats != nullptr);

    vf::grid_winding::SubgridDistanceStats aggregated;

    struct GpuKdAdapter
    {
        GridImp&                           grid;
        uint                               index{ 0 };
        int                                dir{ 0 };
        vf::grid_winding::SubgridDistanceStats* stats{ nullptr };

        bool promoteSolid(const vf::grid_winding::BoundaryLink&) const
        {
            (void)stats;
            return false;
        }

        double existingQ(const vf::grid_winding::BoundaryLink&) const
        {
            if (!grid.hasQIndex(index))
                return -1.0;
            return static_cast<double>(grid.getQValue(index, static_cast<uint>(dir)));
        }

        void setQ(const vf::grid_winding::BoundaryLink&, double q, int patchIndex) const
        {
            if (!grid.hasQIndex(index))
                return;
            grid.setQValue(index, dir, static_cast<real>(q));
            if (patchIndex >= 0)
                grid.setQPatch(index, static_cast<uint>(patchIndex));
        }

        void markSurface(const vf::grid_winding::BoundaryLink&, double, int) const {}

        std::pair<bool, double> windingNumber(const vf::grid_winding::BoundaryLink&) const
        {
            return { false, 0.0 };
        }

        bool handleMissing(const vf::grid_winding::BoundaryLink&,
                           vf::grid_winding::SubgridDistanceStats::MissingLink&) const
        {
            return false;
        }
    };

#pragma omp parallel
    {
#pragma omp single
        VF_LOG_TRACE("Boundary processing OpenMP team size: {}", omp_get_num_threads());

        vf::grid_winding::SubgridDistanceStats localStats;
        auto*                                  statsPtr = collectStats ? &localStats : nullptr;
        GpuKdAdapter                           adapter{ grid, 0U, 0, statsPtr };
        vf::grid_winding::QComputationConfig   cfg;

#pragma omp for schedule(static)
        for (int idx = 0; idx < static_cast<int>(boundaryIndices.size()); ++idx) {
            const uint index = boundaryIndices[static_cast<std::size_t>(idx)];
            if (grid.getFieldEntry(index) != vf::gpu::BC_SOLID)
                continue;

            real ox = 0.0;
            real oy = 0.0;
            real oz = 0.0;
            grid.transIndexToCoords(index, ox, oy, oz);
            vf::grid_winding::Vec3 origin{ static_cast<double>(ox), static_cast<double>(oy),
                                           static_cast<double>(oz) };

            for (int dir = dirStart; dir <= dirEnd; ++dir) {
                const auto& direction = grid.distribution.directions[dir];
                const int   cx        = direction[0];
                const int   cy        = direction[1];
                const int   cz        = direction[2];
                if (cx == 0 && cy == 0 && cz == 0)
                    continue;

                const double neighbourX = origin.X1() + dx * static_cast<double>(cx);
                const double neighbourY = origin.X2() + dx * static_cast<double>(cy);
                const double neighbourZ = origin.X3() + dx * static_cast<double>(cz);

                const uint neighbourIndex = grid.transCoordToIndex(static_cast<real>(neighbourX),
                                                                   static_cast<real>(neighbourY),
                                                                   static_cast<real>(neighbourZ));
                if (neighbourIndex == INVALID_INDEX)
                    continue;

                const char neighbourType = grid.getFieldEntry(neighbourIndex);
                if (!vf::gpu::grid_winding::detail::isSolidNeighbour(neighbourType))
                    continue;

                if (statsPtr)
                    ++statsPtr->totalLinks;

                adapter.index = index;
                adapter.dir   = dir;

                vf::grid_winding::BoundaryLink link;
                link.origin    = origin;
                link.neighbour = vf::grid_winding::Vec3{ neighbourX, neighbourY, neighbourZ };
                link.direction = dir;

                vf::grid_winding::processBoundaryLinkKdTree(kdTrees, adapter, link, statsPtr, cfg);
            }
        }

        if (collectStats) {
#pragma omp critical
            {
                mergeSubgridStats(aggregated, localStats);
            }
        }
    }

    if (collectStats)
        mergeSubgridStats(*stats, aggregated);

    const double boundarySeconds = boundaryTimer.getCurrentRuntimeInSeconds();
    if (collectStats) {
        VF_LOG_TRACE("Boundary link processing finished: hits={}, missing={}, promoted={}, time={} sec",
                     aggregated.hitLinks,
                     aggregated.missingLinks.size(),
                     aggregated.promotedFluidNodes,
                     boundarySeconds);
    } else {
        VF_LOG_TRACE("Boundary link processing finished (no stats): processed {} nodes in {} sec",
                     boundaryIndices.size(),
                     boundarySeconds);
    }
#else
    (void)parallelize;
    (void)parallelThreshold;
    vf::basics::Timer boundaryTimer;
    boundaryTimer.start();

    const double dx       = static_cast<double>(grid.getDelta());
    const int    dirStart = grid.getStartDirection();
    const int    dirEnd   = grid.getEndDirection();

    vf::grid_winding::SubgridDistanceStats localStats;
    auto*                                  statsPtr = stats ? &localStats : nullptr;
    GpuKdAdapter                           adapter{ grid, 0U, 0, statsPtr };
    vf::grid_winding::QComputationConfig   cfg;

    for (const uint index : boundaryIndices) {
        if (grid.getFieldEntry(index) != vf::gpu::BC_SOLID)
            continue;

        real ox = 0.0;
        real oy = 0.0;
        real oz = 0.0;
        grid.transIndexToCoords(index, ox, oy, oz);
        vf::grid_winding::Vec3 origin{ static_cast<double>(ox), static_cast<double>(oy),
                                       static_cast<double>(oz) };

        for (int dir = dirStart; dir <= dirEnd; ++dir) {
            const auto& direction = grid.distribution.directions[dir];
            const int   cx        = direction[0];
            const int   cy        = direction[1];
            const int   cz        = direction[2];
            if (cx == 0 && cy == 0 && cz == 0)
                continue;

            const double neighbourX = origin.X1() + dx * static_cast<double>(cx);
            const double neighbourY = origin.X2() + dx * static_cast<double>(cy);
            const double neighbourZ = origin.X3() + dx * static_cast<double>(cz);

            const uint neighbourIndex = grid.transCoordToIndex(static_cast<real>(neighbourX),
                                                               static_cast<real>(neighbourY),
                                                               static_cast<real>(neighbourZ));
            if (neighbourIndex == INVALID_INDEX)
                continue;

            const char neighbourType = grid.getFieldEntry(neighbourIndex);
            if (!vf::gpu::grid_winding::detail::isSolidNeighbour(neighbourType))
                continue;

            if (statsPtr)
                ++statsPtr->totalLinks;

            adapter.index = index;
            adapter.dir   = dir;

            vf::grid_winding::BoundaryLink link;
            link.origin    = origin;
            link.neighbour = vf::grid_winding::Vec3{ neighbourX, neighbourY, neighbourZ };
            link.direction = dir;

            vf::grid_winding::processBoundaryLinkKdTree(kdTrees, adapter, link, statsPtr, cfg);
        }
    }

    if (stats && statsPtr) {
        mergeSubgridStats(*stats, localStats);
    }

    const double boundarySeconds = boundaryTimer.getCurrentRuntimeInSeconds();
    if (stats) {
        VF_LOG_TRACE("Boundary link processing finished: hits={}, missing={}, promoted={}, time={} sec",
                     stats->hitLinks,
                     stats->missingLinks.size(),
                     stats->promotedFluidNodes,
                     boundarySeconds);
    } else {
        VF_LOG_TRACE("Boundary link processing finished (serial) in {} sec", boundarySeconds);
    }
#endif
}

} // namespace vf::grid_winding::gpu
