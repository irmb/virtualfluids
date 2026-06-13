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
//! \addtogroup gpu_io io
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Hussein Alihussein
//! \brief Implements GPU grid-winding diagnostics collection and VTK output helpers.
//! \note Generated with the assistance of OpenAI Codex (GPT-5). Reviewed and adapted by author.
//=======================================================================================
#include "GridWindingDiagnosticsGpu.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <memory>
#include <utility>

#include <basics/container/CbArray3D.h>
#include <basics/Timer/Timer.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <basics/geometry3d/winding/GridWindingWriting.h>
#include <basics/geometry3d/winding/GridWindingFastDefaults.h>
#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>

#include <logger/Logger.h>

#include <parallel/Communicator.h>

#include "grid/Grid.h"
#include "grid/GridImp.h"
#include "grid/NodeValues.h"

#if defined(VF_HAS_FAST_WINDING)
#include "utilities/GridWindingBoundaryGpuCore.h"

#include <igl/FastWindingNumberForSoups.h>
#include <igl/parallel_for.h>
#endif

using vf::grid_winding::QLineCollection;
namespace vf::gpu {
void GridWindingWriter::writeMissingLinks(const vf::grid_winding::MissingLinkCollection& collection,
    const std::string&                              basePath)
    {
        vf::grid_winding::writeMissingLinks(collection, basePath);
    }
    
void GridWindingWriter::writeQLines(const vf::grid_winding::QLineCollection& collection,
    const std::string&                        basePath)
    {
        vf::grid_winding::writeQLines(collection, basePath);
    }

namespace grid_winding
{

void collectQLines(const std::vector<SPtr<Grid>>&          grids,
                   const SPtr<vf::parallel::Communicator>& comm,
                   QLineCollection&                        collection)
{
    collection = QLineCollection{};
    if (grids.empty())
        return;

    std::vector<double> coordsSend;
    std::vector<int>    dirSend;
    std::vector<float>  qSend;
    for (const auto& gridBase : grids) {
        auto gridImp = std::dynamic_pointer_cast<GridImp>(gridBase);
        if (!gridImp)
            continue;

        const double dx       = static_cast<double>(gridImp->getDelta());
        const int    dirStart = gridImp->distribution.dir_start;
        const int    dirEnd   = gridImp->distribution.dir_end;

        const uint nodeCount = gridImp->getSize();
        for (uint index = 0; index < nodeCount; ++index) {
            if (!gridImp->hasQIndex(index))
                continue;

            if (gridImp->getFieldEntry(index) != vf::gpu::BC_SOLID)
                continue;

            real ox = 0.0;
            real oy = 0.0;
            real oz = 0.0;
            gridImp->transIndexToCoords(index, ox, oy, oz);

            for (int dir = dirStart; dir <= dirEnd; ++dir) {
                const real   qValue = gridImp->getQValue(index, static_cast<uint>(dir));
                const double rawQ   = static_cast<double>(qValue);

                if (ub_math::lessEqual(qValue, static_cast<real>(0.0)))
                    continue;

                const double clampedQ = std::clamp(rawQ, 0.0, 1.0);
                if (ub_math::zero(clampedQ))
                    continue;

                const auto& direction = gridImp->distribution.directions[dir];
                const int   cx        = direction[0];
                const int   cy        = direction[1];
                const int   cz        = direction[2];
                if (cx == 0 && cy == 0 && cz == 0)
                    continue;

                const double hx = static_cast<double>(ox) + clampedQ * dx * static_cast<double>(cx);
                const double hy = static_cast<double>(oy) + clampedQ * dx * static_cast<double>(cy);
                const double hz = static_cast<double>(oz) + clampedQ * dx * static_cast<double>(cz);

                coordsSend.push_back(static_cast<double>(ox));
                coordsSend.push_back(static_cast<double>(oy));
                coordsSend.push_back(static_cast<double>(oz));
                coordsSend.push_back(hx);
                coordsSend.push_back(hy);
                coordsSend.push_back(hz);

                dirSend.push_back(dir);
                qSend.push_back(static_cast<float>(clampedQ));
            }
        }
    }

    std::vector<double> coordsRecv;
    std::vector<int>    dirRecv;
    std::vector<float>  qRecv;
    if (comm) {
        comm->allGather(coordsSend, coordsRecv);
        comm->allGather(dirSend, dirRecv);
        comm->allGather(qSend, qRecv);
    } else {
        coordsRecv = std::move(coordsSend);
        dirRecv    = std::move(dirSend);
        qRecv      = std::move(qSend);
    }

    if (comm && !comm->isRoot())
        return;

    const bool hasLineData =
        !qRecv.empty() && coordsRecv.size() >= qRecv.size() * 6;

    if (!hasLineData)
        return;

    collection = vf::grid_winding::buildQLineCollection(
        coordsRecv, dirRecv, qRecv);
}

void runBoundaryDiagnosticsForSurfaces(const std::vector<SPtr<Grid>>&          grids,
                                       const SPtr<GbTriFaceMesh3D>&            buildingSurface,
                                       const SPtr<GbTriFaceMesh3D>&            bodenSurface,
                                       bool                                    writeMissingLinks,
                                       bool                                    writeQLinesFlag,
                                       const std::string&                      path,
                                       const SPtr<vf::parallel::Communicator>& comm,
                                       bool                                    isRoot)
{
#if defined(VF_HAS_FAST_WINDING)
    if (grids.empty())
        return;

    std::vector<GbTriFaceMesh3D*> surfaces;
    surfaces.reserve(2);
    if (buildingSurface && buildingSurface->getTriangles() && !buildingSurface->getTriangles()->empty())
        surfaces.push_back(buildingSurface.get());
    if (bodenSurface && bodenSurface->getTriangles() && !bodenSurface->getTriangles()->empty())
        surfaces.push_back(bodenSurface.get());

    const auto combineSurfaces = [](const SPtr<GbTriFaceMesh3D>& first,
                                    const SPtr<GbTriFaceMesh3D>& second) -> SPtr<GbTriFaceMesh3D> {
        if (!first && !second)
            return nullptr;
        if (!second)
            return first;
        if (!first)
            return second;

        std::unique_ptr<GbTriFaceMesh3D> combined(first->clone());
        auto* combinedNodes = combined->getNodes();
        auto* combinedTris = combined->getTriangles();
        const auto* secondNodes = second->getNodes();
        const auto* secondTris = second->getTriangles();
        if (!combinedNodes || !combinedTris || !secondNodes || !secondTris)
            return first;

        const auto baseIndex = combinedNodes->size();
        combinedNodes->reserve(baseIndex + secondNodes->size());
        for (const auto& node : *secondNodes)
            combinedNodes->push_back(node);

        combinedTris->reserve(combinedTris->size() + secondTris->size());
        for (const auto& tri : *secondTris) {
            combinedTris->emplace_back(
                static_cast<int>(tri.getIndexVertex1() + baseIndex),
                static_cast<int>(tri.getIndexVertex2() + baseIndex),
                static_cast<int>(tri.getIndexVertex3() + baseIndex));
        }

        combined->calculateValues();
        return SPtr<GbTriFaceMesh3D>(combined.release());
    };

    const auto isFluidType = [](char type) {
        switch (type) {
        case FLUID:
        case FLUID_CFC:
        case FLUID_CFF:
        case FLUID_FCC:
        case FLUID_FCF:
            return true;
        default:
            return false;
        }
    };

    const auto isBoundaryType = [](char type) {
        switch (type) {
        case BC_SOLID:
        case BC_VELOCITY:
        case BC_PRESSURE:
        case BC_SLIP:
        case BC_STRESS:
        case BC_OUTFLOW:
            return true;
        default:
            return false;
        }
    };

    const auto combinedSurface = combineSurfaces(buildingSurface, bodenSurface);
    if (isRoot) {
        const auto nodeFlagProvider = [isFluidType, isBoundaryType](char type) {
            return visualization::NodeFlags{
                isFluidType(type) ? 1.0 : 0.0,
                isBoundaryType(type) ? 1.0 : 0.0
            };
        };
        visualization::exportGridVisualization(
            grids, combinedSurface, path + "/vtk", nodeFlagProvider);
    }

    std::vector<GbTriFaceMesh3D*> kdMeshes;
    kdMeshes.reserve(2);
    if (buildingSurface)
        kdMeshes.push_back(buildingSurface.get());
    if (bodenSurface)
        kdMeshes.push_back(bodenSurface.get());

    runBoundaryDiagnostics(grids,
                           surfaces,
                           kdMeshes,
                           writeMissingLinks,
                           writeQLinesFlag,
                           path,
                           comm,
                           isRoot);
#else
    (void)buildingSurface;
    (void)bodenSurface;
    (void)writeMissingLinks;
    (void)writeQLinesFlag;
    (void)path;
    (void)comm;

    if (!grids.empty() && isRoot) {
        VF_LOG_INFO("Grid-winding diagnostics disabled (VF_HAS_FAST_WINDING not enabled).");
    }
#endif
}

} // namespace vf::gpu::grid_winding

#if defined(VF_HAS_FAST_WINDING)
namespace
{
using NodeFlags = vf::gpu::grid_winding::visualization::NodeFlags;
using NodeFlagProvider = vf::gpu::grid_winding::visualization::NodeFlagProvider;

NodeFlags evaluateNodeFlags(char type, const NodeFlagProvider& provider)
{
    if (provider)
        return provider(type);
    return NodeFlags{};
}

std::vector<float> computeFastWindingField(GridImp& grid,
                                           const SPtr<GbTriFaceMesh3D>& surface,
                                           float accuracyScale = vf::grid_winding::FastWindingDefaultAccuracyScale)
{
    std::vector<float> winding(grid.getSize(), std::numeric_limits<float>::quiet_NaN());

    if (!surface)
        return winding;

    const auto* nodes = surface->getNodes();
    const auto* tris  = surface->getTriangles();
    if (!nodes || !tris || nodes->empty() || tris->empty())
        return winding;

    using FastVector     = igl::FastWindingNumber::HDK_Sample::UT_Vector3T<float>;
    using FastSolidAngle = igl::FastWindingNumber::HDK_Sample::UT_SolidAngle<float, float>;

    std::vector<FastVector> meshPositions;
    meshPositions.reserve(nodes->size());
    for (const auto& node : *nodes) {
        FastVector p{};
        p[0] = static_cast<float>(node.x);
        p[1] = static_cast<float>(node.y);
        p[2] = static_cast<float>(node.z);
        meshPositions.push_back(p);
    }

    std::vector<int> meshIndices;
    meshIndices.reserve(tris->size() * 3);
    for (const auto& tri : *tris) {
        meshIndices.push_back(static_cast<int>(tri.getIndexVertex1()));
        meshIndices.push_back(static_cast<int>(tri.getIndexVertex2()));
        meshIndices.push_back(static_cast<int>(tri.getIndexVertex3()));
    }

    //! libigl fast-winding implementation (UT_SolidAngle) follows:
    //! \ref <a href="https://doi.org/10.1145/3197517.3201337"><b>[ G. Barill et al. (2018), DOI:10.1145/3197517.3201337 ]</b></a>
    FastSolidAngle fastWinding;
    fastWinding.init(static_cast<int>(tris->size()), meshIndices.data(), static_cast<int>(meshPositions.size()),
                     meshPositions.data(), vf::grid_winding::FastWindingDefaultExpansionOrder);

    const float invFourPi = 1.0f / (4.0f * static_cast<float>(std::acos(-1.0)));

    const uint totalSamples = grid.getSize();
    if (totalSamples > 0) {
        igl::parallel_for(static_cast<int>(totalSamples), [&](int index) {
            real x{ 0.0f };
            real y{ 0.0f };
            real z{ 0.0f };
            grid.transIndexToCoords(static_cast<uint>(index), x, y, z);

            FastVector sample{};
            sample[0] = static_cast<float>(x);
            sample[1] = static_cast<float>(y);
            sample[2] = static_cast<float>(z);

            // `accuracyScale` is the libigl fast-winding accuracy/speed parameter from the cited method.
            const float angle        = fastWinding.computeSolidAngle(sample, accuracyScale);
            const float windingValue = angle * invFourPi;
            winding[static_cast<std::size_t>(index)] = windingValue;
        }, static_cast<size_t>(1024));
    }

    return winding;
}

void writeGridVisualization(GridImp& grid, const std::vector<float>& winding, const std::string& baseName,
                            NodeFlagProvider flagProvider)
{
    using ::makeUbTuple;
    const uint nx = grid.getNumberOfNodesX();
    const uint ny = grid.getNumberOfNodesY();
    const uint nz = grid.getNumberOfNodesZ();

    const uint64_t planeSize  = static_cast<uint64_t>(nx) * static_cast<uint64_t>(ny);
    const uint     chunkSize  = 20000000;
    const uint     chunkSizeZ = std::max<uint>(1, static_cast<uint>(chunkSize / std::max<uint64_t>(1, planeSize)));

    for (uint startZ = 0, endZ = std::min<uint>(nz, chunkSizeZ), part = 0; startZ < nz;
         startZ = endZ, endZ = std::min<uint>(nz, endZ + chunkSizeZ), ++part) {

        std::vector<UbTupleFloat3>       nodes;
        std::vector<UbTupleUInt8>        cells;
        std::vector<std::string>         nodedatanames{ "field_type", "is_fluid", "is_boundary", "winding_number",
                                                "sparse_id",  "matrix_id" };
        std::vector<std::vector<double>> nodedata(nodedatanames.size());

        CbArray3D<int> nodeNumbers(nx, ny, nz, -1);
        int            currentNodeIndex = 0;

        for (uint xIndex = 0; xIndex < nx; ++xIndex) {
            for (uint yIndex = 0; yIndex < ny; ++yIndex) {
                for (uint zIndex = startZ; zIndex < endZ; ++zIndex) {
                    const uint matrixIndex = nx * ny * zIndex + nx * yIndex + xIndex;

                    real x{ 0.0f };
                    real y{ 0.0f };
                    real z{ 0.0f };
                    grid.transIndexToCoords(matrixIndex, x, y, z);

                    nodeNumbers(xIndex, yIndex, zIndex) = currentNodeIndex++;
                    nodes.emplace_back(UbTupleFloat3(float(x), float(y), float(z)));

                    const char      type  = grid.getFieldEntry(matrixIndex);
                    nodedata[0].push_back(static_cast<double>(static_cast<int>(type)));
                    const NodeFlags flags = evaluateNodeFlags(type, flagProvider);
                    nodedata[1].push_back(flags.isFluid);
                    nodedata[2].push_back(flags.isBoundary);
                    double windingValue = std::numeric_limits<double>::quiet_NaN();
                    if (!winding.empty() && matrixIndex < winding.size())
                        windingValue = static_cast<double>(winding[matrixIndex]);
                    nodedata[3].push_back(windingValue);
                    nodedata[4].push_back(static_cast<double>(grid.getSparseIndex(matrixIndex)));
                    nodedata[5].push_back(static_cast<double>(matrixIndex));
                }
            }
        }

        for (uint xIndex = 0; xIndex + 1 < nx; ++xIndex) {
            for (uint yIndex = 0; yIndex + 1 < ny; ++yIndex) {
                for (uint zIndex = startZ; zIndex + 1 < endZ; ++zIndex) {
                    const int SWB = nodeNumbers(xIndex, yIndex, zIndex);
                    const int SEB = nodeNumbers(xIndex + 1, yIndex, zIndex);
                    const int NEB = nodeNumbers(xIndex + 1, yIndex + 1, zIndex);
                    const int NWB = nodeNumbers(xIndex, yIndex + 1, zIndex);
                    const int SWT = nodeNumbers(xIndex, yIndex, zIndex + 1);
                    const int SET = nodeNumbers(xIndex + 1, yIndex, zIndex + 1);
                    const int NET = nodeNumbers(xIndex + 1, yIndex + 1, zIndex + 1);
                    const int NWT = nodeNumbers(xIndex, yIndex + 1, zIndex + 1);

                    if (SWB < 0 || SEB < 0 || NEB < 0 || NWB < 0 || SWT < 0 || SET < 0 || NET < 0 || NWT < 0)
                        continue;

                    cells.emplace_back(makeUbTuple(uint(SWB), uint(SEB), uint(NEB), uint(NWB), uint(SWT), uint(SET),
                                                   uint(NET), uint(NWT)));
                }
            }
        }

        const std::string filename = baseName + "_part_" + std::to_string(part);
        VF_LOG_INFO("Writing grid visualization to {}", filename);
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(filename, nodes, cells, nodedatanames, nodedata);
    }
}

} // namespace

namespace grid_winding
{
void runBoundaryDiagnostics(const std::vector<SPtr<Grid>>&                          grids,
                            const std::vector<GbTriFaceMesh3D*>&                    surfaces,
                            std::vector<GbTriFaceMesh3D*>&                          kdMeshes,
                            bool                                                    writeMissingLinks,
                            bool                                                    writeQLinesFlag,
                            const std::string&                                      path,
                            const SPtr<parallel::Communicator>&                 comm,
                            bool                                                    isRoot)
{
    if (grids.empty())
        return;

    vf::gpu::grid_winding::BoundaryProcessingConfig boundaryConfig;

    vf::grid_winding::SubgridDistanceStats subgridStats;

    if (!surfaces.empty()) {
        vf::basics::Timer statsTimer;
        statsTimer.start();

        for (const auto& gridBase : grids) {
            auto gridImp = std::dynamic_pointer_cast<GridImp>(gridBase);
            if (!gridImp)
                continue;

            if (boundaryConfig.rebuildBoundaryQIndices)
                gridImp->rebuildBoundaryQIndices();

            if (boundaryConfig.ensureQStorage)
                gridImp->ensureQStorageAllocated();

            vf::grid_winding::SubgridDistanceStats localStats;
            vf::gpu::grid_winding::processSolidBoundaryLinksMultithreaded(
                *gridImp, surfaces, &localStats, boundaryConfig.parallelize, boundaryConfig.parallelThreshold,
                &kdMeshes);
            vf::gpu::grid_winding::mergeSubgridStats(subgridStats, localStats);

            if (boundaryConfig.fillMissingQs)
                gridImp->fillMissingQsAlongSolidNeighbours(boundaryConfig.missingQValue);
        }

        VF_LOG_INFO("Subgrid statistics collected in {} sec", statsTimer.getCurrentRuntimeInSeconds());
    } else {
        VF_LOG_INFO("Subgrid statistics skipped (no surface data available).");
    }

    if (writeMissingLinks) {
        vf::basics::Timer missingLinksTimer;
        missingLinksTimer.start();
        const auto missingLinks = vf::grid_winding::collectMissingLinks(subgridStats, comm);
        if (isRoot) {
            VF_LOG_INFO("Missing link aggregation completed in {} sec",
                        missingLinksTimer.getCurrentRuntimeInSeconds());
            std::error_code bcDirStatus;
            std::filesystem::create_directories(path + "/bc", bcDirStatus);
            vf::basics::Timer missingLinksWriteTimer;
            missingLinksWriteTimer.start();
            GridWindingWriter::writeMissingLinks(missingLinks, path + "/bc/missing_boundary_links");
            VF_LOG_INFO("Missing link file written in {} sec",
                        missingLinksWriteTimer.getCurrentRuntimeInSeconds());
        }
    } else if (isRoot) {
        VF_LOG_INFO("Missing link export disabled; skipping file generation.");
    }

    if (writeQLinesFlag) {
        vf::basics::Timer qCollectTimer;
        qCollectTimer.start();
        vf::grid_winding::QLineCollection qCollection;
        vf::gpu::grid_winding::collectQLines(grids, comm, qCollection);
        if (isRoot) {
            const bool hasLines = !qCollection.empty();
            if (hasLines) {
                VF_LOG_INFO("Q-line collection completed in {} sec",
                            qCollectTimer.getCurrentRuntimeInSeconds());
                vf::basics::Timer qWriteTimer;
                qWriteTimer.start();
                GridWindingWriter::writeQLines(qCollection, path + "/bc/boundary_qs");
                VF_LOG_INFO("Q-line file written in {} sec", qWriteTimer.getCurrentRuntimeInSeconds());
            } else {
                VF_LOG_INFO("Q-line export disabled; skipping file generation.");
            }
        }
    } else if (isRoot) {
        VF_LOG_INFO("Q-line export disabled; skipping collection.");
    }
}

namespace visualization
{
void exportGridVisualization(const std::vector<SPtr<Grid>>& grids,
                             const SPtr<GbTriFaceMesh3D>&   surface,
                             const std::string&            outputPath,
                             NodeFlagProvider              flagProvider)
{
    if (grids.empty())
        return;

    std::error_code ec;
    std::filesystem::create_directories(outputPath, ec);
    if (ec)
        VF_LOG_WARNING("Unable to create output directory {}: {}", outputPath, ec.message());

    for (std::size_t level = 0; level < grids.size(); ++level) {
        auto gridImp = std::dynamic_pointer_cast<GridImp>(grids[level]);
        if (!gridImp)
            continue;

        std::vector<float> windingValues;
        if (surface && surface->getTriangles() && !surface->getTriangles()->empty()) {
            VF_LOG_TRACE("Computing winding numbers for grid level {}", level);
            vf::basics::Timer windingTimer;
            windingTimer.start();
            windingValues = computeFastWindingField(*gridImp, surface);
            VF_LOG_INFO("Grid level {} winding number computed in {} sec", level,
                        windingTimer.getCurrentRuntimeInSeconds());
        } else {
            VF_LOG_INFO("Grid level {} winding computation skipped (no surface data)", level);
        }

        const std::string fileBase = outputPath + "/grid_level_" + std::to_string(level);
        VF_LOG_TRACE("Writing grid visualization for level {}", level);
        vf::basics::Timer writingTimer;
        writingTimer.start();
        writeGridVisualization(*gridImp, windingValues, fileBase, flagProvider);
        VF_LOG_INFO("Grid level {} visualization written in {} sec", level,
                    writingTimer.getCurrentRuntimeInSeconds());
    }
}
} // namespace visualization
} // namespace vf::gpu::grid_winding
}
#endif

//! \}
