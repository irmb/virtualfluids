//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /    \     |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /      \    |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\    \   |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____    \  |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \____\ |_______|
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
//! \author Hussein Alihussein
//! \note CPU grid-winding diagnostics observer, introduced with assistance from GPT-5.1.
//=======================================================================================

#include "GridWindingDiagnosticsSimulationObserver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <basics/Timer/Timer.h>
#include <logger/Logger.h>

#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>
#include <basics/geometry3d/winding/GridWindingFastDefaults.h>
#include <basics/geometry3d/winding/GridWindingWriting.h>

#include <cpu/core/Interactors/D3Q27GridWindingInteractor.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>

#include <parallel/Communicator.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include "BCArray3D.h"
#include "Block3D.h"
#include "CbArray3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "BCSet.h"
#include "UbScheduler.h"
#include "WbWriter.h"

#include <filesystem>

#if defined(VF_HAS_FAST_WINDING)
#include <igl/FastWindingNumberForSoups.h>
#endif

namespace
{
SPtr<GbTriFaceMesh3D> combineSurfaces(const std::vector<SPtr<GbTriFaceMesh3D>> &surfaces)
{
    std::vector<SPtr<GbTriFaceMesh3D>> filtered;
    filtered.reserve(surfaces.size());
    for (const auto &s : surfaces) {
        if (!s)
            continue;
        const auto *nodes = s->getNodes();
        const auto *tris  = s->getTriangles();
        if (!nodes || !tris || nodes->empty() || tris->empty())
            continue;
        filtered.push_back(s);
    }

    if (filtered.empty())
        return nullptr;
    if (filtered.size() == 1)
        return filtered.front();

    auto combined  = std::make_shared<GbTriFaceMesh3D>();
    auto *allNodes = combined->getNodes();
    auto *allTris  = combined->getTriangles();

    for (const auto &surface : filtered) {
        const auto *srcNodes = surface->getNodes();
        const auto *srcTris  = surface->getTriangles();
        if (!srcNodes || !srcTris)
            continue;

        const int nodeOffset = static_cast<int>(allNodes->size());
        allNodes->reserve(allNodes->size() + srcNodes->size());
        allTris->reserve(allTris->size() + srcTris->size());

        allNodes->insert(allNodes->end(), srcNodes->begin(), srcNodes->end());
        for (const auto &tri : *srcTris) {
            allTris->emplace_back(tri.getIndexVertex1() + nodeOffset,
                                  tri.getIndexVertex2() + nodeOffset,
                                  tri.getIndexVertex3() + nodeOffset);
        }
    }

    combined->calculateValues();
    return combined;
}

#if defined(VF_HAS_FAST_WINDING)
double gridCellSizeHint(const SPtr<Grid3D> &grid)
{
    if (!grid)
        return 0.0;

    const int coarsestLevel = grid->getCoarsestInitializedLevel();
    std::vector<SPtr<Block3D>> blocks;
    grid->getBlocks(coarsestLevel, blocks);
    if (blocks.empty())
        return 0.0;

    return static_cast<double>(grid->getDeltaX(blocks.front()));
}

struct FastWindingEvaluator
{
    using FastVector     = igl::FastWindingNumber::HDK_Sample::UT_Vector3T<float>;
    using FastSolidAngle = igl::FastWindingNumber::HDK_Sample::UT_SolidAngle<float, float>;

    std::array<double, 3> meshMin{};
    std::array<double, 3> meshMax{};
    std::vector<FastVector> meshPositions{};
    std::vector<int> meshIndices{};
    FastSolidAngle         fastWinding{};
    float                  accuracyScale{ vf::grid_winding::FastWindingDefaultAccuracyScale };
    float                  invFourPi{ 1.0f };
    float                  clampMax{ 1.0f };
    double                 faceEpsBase{ 0.0 };
    bool                   ready{ false };

    [[nodiscard]] float computeClamped(double x, double y, double z) const
    {
        if (!ready)
            return std::numeric_limits<float>::quiet_NaN();

        FastVector query{};
        query[0] = static_cast<float>(x);
        query[1] = static_cast<float>(y);
        query[2] = static_cast<float>(z);

        // `accuracyScale` is the libigl fast-winding accuracy/speed parameter from the cited method.
        const float winding = fastWinding.computeSolidAngle(query, accuracyScale) * invFourPi;
        return std::clamp(winding, 0.0f, clampMax);
    }

    [[nodiscard]] bool computeClampedIfInside(double x, double y, double z, double faceEps, float &out) const
    {
        if (!ready)
            return false;

        if (ub_math::less(x, meshMin[0] - faceEps) || ub_math::greater(x, meshMax[0] + faceEps) ||
            ub_math::less(y, meshMin[1] - faceEps) || ub_math::greater(y, meshMax[1] + faceEps) ||
            ub_math::less(z, meshMin[2] - faceEps) || ub_math::greater(z, meshMax[2] + faceEps))
            return false;

        out = computeClamped(x, y, z);
        return true;
    }
};

bool initFastWindingEvaluator(FastWindingEvaluator &eval, const SPtr<GbTriFaceMesh3D> &surface,
                              const SPtr<Grid3D> &grid,
                              float accuracyScale = vf::grid_winding::FastWindingDefaultAccuracyScale,
                              float tolerance = vf::grid_winding::FastWindingDefaultTolerance)
{
    if (!surface)
        return false;

    const auto *nodes = surface->getNodes();
    const auto *tris  = surface->getTriangles();
    if (!nodes || !tris || nodes->empty() || tris->empty())
        return false;

    eval.meshPositions.clear();
    eval.meshPositions.reserve(nodes->size());
    for (const auto &node : *nodes) {
        FastWindingEvaluator::FastVector p{};
        p[0] = static_cast<float>(node.x);
        p[1] = static_cast<float>(node.y);
        p[2] = static_cast<float>(node.z);
        eval.meshPositions.push_back(p);
    }

    eval.meshIndices.clear();
    eval.meshIndices.reserve(tris->size() * 3);
    for (const auto &tri : *tris) {
        eval.meshIndices.push_back(static_cast<int>(tri.getIndexVertex1()));
        eval.meshIndices.push_back(static_cast<int>(tri.getIndexVertex2()));
        eval.meshIndices.push_back(static_cast<int>(tri.getIndexVertex3()));
    }

    //! libigl fast-winding implementation (UT_SolidAngle) follows:
    //! \ref <a href="https://doi.org/10.1145/3197517.3201337"><b>[ G. Barill et al. (2018), DOI:10.1145/3197517.3201337 ]</b></a>
    eval.fastWinding.init(static_cast<int>(tris->size()), eval.meshIndices.data(),
                          static_cast<int>(eval.meshPositions.size()), eval.meshPositions.data(),
                          vf::grid_winding::FastWindingDefaultExpansionOrder);

    const float effectiveTolerance = ub_math::zero(tolerance) ? vf::grid_winding::FastWindingDefaultTolerance
                                                              : tolerance;

    eval.meshMin = { surface->getX1Minimum(), surface->getX2Minimum(), surface->getX3Minimum() };
    eval.meshMax = { surface->getX1Maximum(), surface->getX2Maximum(), surface->getX3Maximum() };
    eval.accuracyScale = accuracyScale;
    eval.invFourPi     = 1.0f / (4.0f * static_cast<float>(std::acos(-1.0)));
    eval.clampMax      = 1.0f + 10.0f * effectiveTolerance;
    eval.faceEpsBase   = vf::grid_winding::computeFaceEpsilon(eval.meshMin, eval.meshMax, gridCellSizeHint(grid));
    eval.ready         = true;
    return true;
}

void populateMissingLinkWinding(vf::grid_winding::MissingLinkCollection &collection,
                                const FastWindingEvaluator &eval)
{
    if (collection.total == 0 || collection.nodes.size() < collection.total * 2)
        return;

    collection.windingNumbers.resize(collection.total, std::numeric_limits<float>::quiet_NaN());
    for (std::size_t idx = 0; idx < collection.total; ++idx) {
        const auto &node = collection.nodes[idx * 2 + 1];
        const float x    = val<1>(node);
        const float y    = val<2>(node);
        const float z    = val<3>(node);
        collection.windingNumbers[idx] = eval.computeClamped(static_cast<double>(x),
                                                             static_cast<double>(y),
                                                             static_cast<double>(z));
    }
}
#endif
} // namespace

namespace vf::grid_winding::cpu
{
void finalizeMissingLinkCollection(vf::grid_winding::MissingLinkCollection &collection,
                                   const SPtr<GbTriFaceMesh3D>            &surface,
                                   const SPtr<Grid3D>                     &grid)
{
    if (collection.empty())
        return;

#if defined(VF_HAS_FAST_WINDING)
    if (surface && grid) {
        FastWindingEvaluator eval;
        if (initFastWindingEvaluator(eval, surface, grid))
            populateMissingLinkWinding(collection, eval);
    }
#endif
    if (surface)
        vf::grid_winding::populateMissingLinkEndpointKind(collection, surface);

    if (!collection.endpointKind.empty() && collection.endpointKind.size() == collection.total) {
        vf::grid_winding::MissingLinkCollection filtered;
        filtered.nodes.reserve(collection.nodes.size());
        filtered.lines.reserve(collection.lines.size());
        filtered.directions.reserve(collection.directions.size());
        filtered.windingNumbers.reserve(collection.windingNumbers.size());
        filtered.endpointKind.reserve(collection.endpointKind.size());

        for (std::size_t idx = 0; idx < collection.total; ++idx) {
            if (!ub_math::zero(collection.endpointKind[idx]))
                continue;

            const int nodeIndex = static_cast<int>(filtered.nodes.size());
            filtered.nodes.emplace_back(collection.nodes[idx * 2 + 0]);
            filtered.nodes.emplace_back(collection.nodes[idx * 2 + 1]);
            filtered.lines.emplace_back(nodeIndex, nodeIndex + 1);

            if (idx < collection.directions.size())
                filtered.directions.push_back(collection.directions[idx]);
            if (idx < collection.windingNumbers.size())
                filtered.windingNumbers.push_back(collection.windingNumbers[idx]);
            filtered.endpointKind.push_back(0.0f);
        }

        filtered.total = filtered.lines.size();
        collection = std::move(filtered);
    }
}
} // namespace vf::grid_winding::cpu

GridWindingDiagnosticsSimulationObserver::GridWindingDiagnosticsSimulationObserver(
    SPtr<Grid3D>                                    grid,
    SPtr<UbScheduler>                               scheduler,
    const std::vector<SPtr<GbTriFaceMesh3D>>       &surfaces,
    const SPtr<BCSet>                              &bcSet,
    const std::vector<SPtr<BC>>                    &boundaryAdapters,
    const SPtr<D3Q27GridWindingInteractor>         &interactor,
    const std::string                              &path,
    const SPtr<vf::parallel::Communicator>         &comm,
    const vf::grid_winding::cpu::BoundaryProcessingConfig &config)
    : SimulationObserver(std::move(grid), std::move(scheduler))
    , path_(path)
    , comm_(comm)
    , surface_(combineSurfaces(surfaces))
    , bcSet_(bcSet)
    , boundaryAdapters_(boundaryAdapters)
    , interactor_(interactor)
    , config_(config)
{
}

GridWindingDiagnosticsSimulationObserver::GridWindingDiagnosticsSimulationObserver(
    SPtr<Grid3D>                                    grid,
    SPtr<UbScheduler>                               scheduler,
    const SPtr<GbTriFaceMesh3D>                    &surface,
    const SPtr<BCSet>                              &bcSet,
    const std::vector<SPtr<BC>>                    &boundaryAdapters,
    const SPtr<D3Q27GridWindingInteractor>         &interactor,
    const std::string                              &path,
    const SPtr<vf::parallel::Communicator>         &comm,
    const vf::grid_winding::cpu::BoundaryProcessingConfig &config)
    : GridWindingDiagnosticsSimulationObserver(
          std::move(grid),
          std::move(scheduler),
          std::vector<SPtr<GbTriFaceMesh3D>>{ surface },
          bcSet,
          boundaryAdapters,
          interactor,
          path,
          comm,
          config)
{
}

void GridWindingDiagnosticsSimulationObserver::update(real step)
{
    if (!scheduler || !scheduler->isDue(step))
        return;

    if (diagnosticsComputed_)
        return;

    diagnosticsComputed_ = true;

    if (!grid || !surface_ || !bcSet_) {
        VF_LOG_WARNING("GridWindingDiagnosticsSimulationObserver::update - missing grid/geometry/BCSet.");
        return;
    }

    auto comm   = comm_ ? comm_ : vf::parallel::Communicator::getInstance();
    const bool isRoot = !comm || comm->isRoot();

    vf::grid_winding::cpu::BoundaryProcessingConfig effectiveConfig = config_;

    vf::basics::Timer timer;
    timer.start();

    if (interactor_) {
        interactor_->computeBoundaryData(surface_, grid, bcSet_, boundaryAdapters_, static_cast<double>(step),
                                         comm, path_, effectiveConfig);
        const auto &result = interactor_->getBoundaryResult();
        stats_             = result.stats;
        missingLinks_      = result.missingLinks;
        qLines_            = result.qLines;
    } else {
        auto result = vf::grid_winding::cpu::processBoundaryData(
            surface_, grid, bcSet_, boundaryAdapters_, nullptr, static_cast<double>(step), comm, path_,
            effectiveConfig);
        stats_        = result.stats;
        missingLinks_ = std::move(result.missingLinks);
        qLines_       = std::move(result.qLines);
    }

    if (isRoot) {
        VF_LOG_INFO("GridWindingDiagnosticsSimulationObserver - diagnostics completed in {} sec",
                    timer.getCurrentRuntimeInSeconds());
    }

    ensureOutputDirectory();
    if (config_.writeMissingLinks)
        writeMissingLinkReport();
    if (config_.writeQLines)
        writeQReports();
    writeWindingVolumeReport(step);
}

void GridWindingDiagnosticsSimulationObserver::ensureOutputDirectory() const
{
    if (path_.empty())
        return;

    std::error_code ec;
    std::filesystem::create_directories(path_ + "/bc", ec);
    if (ec) {
        VF_LOG_WARNING("GridWindingDiagnosticsSimulationObserver - unable to create directory {}: {}",
                       path_ + "/bc", ec.message());
    }
}

void GridWindingDiagnosticsSimulationObserver::writeMissingLinkReport() const
{
    auto comm = comm_ ? comm_ : vf::parallel::Communicator::getInstance();
    const bool isRoot = !comm || comm->isRoot();

    if (!isRoot || missingLinks_.empty())
        return;

    auto collection = missingLinks_;
    vf::grid_winding::cpu::finalizeMissingLinkCollection(collection, surface_, grid);
    if (collection.empty())
        return;

    vf::basics::Timer writeTimer;
    writeTimer.start();
    vf::grid_winding::writeMissingLinks(collection, path_ + "/bc/missing_boundary_links");
    VF_LOG_INFO("GridWindingDiagnosticsSimulationObserver - missing-link file written in {} sec",
                writeTimer.getCurrentRuntimeInSeconds());
}

void GridWindingDiagnosticsSimulationObserver::writeQReports() const
{
    auto comm = comm_ ? comm_ : vf::parallel::Communicator::getInstance();
    const bool isRoot = !comm || comm->isRoot();

    if (!isRoot)
        return;

    if (qLines_.empty()) {
        VF_LOG_INFO("GridWindingDiagnosticsSimulationObserver - no Q-lines to export.");
        return;
    }

    vf::basics::Timer writeTimer;
    writeTimer.start();
    vf::grid_winding::writeQLines(qLines_, path_ + "/bc/boundary_qs");
    VF_LOG_INFO("GridWindingDiagnosticsSimulationObserver - Q-line file written in {} sec",
                writeTimer.getCurrentRuntimeInSeconds());
}

void GridWindingDiagnosticsSimulationObserver::writeWindingVolumeReport(real step) const
{
    auto comm = comm_ ? comm_ : vf::parallel::Communicator::getInstance();
    if (!grid || !surface_)
        return;

#if defined(VF_HAS_FAST_WINDING)
    FastWindingEvaluator eval;
    const bool hasWinding = initFastWindingEvaluator(eval, surface_, grid);
#endif

    const int gridRank     = comm ? comm->getProcessID() : 0;
    const int minInitLevel = grid->getCoarsestInitializedLevel();
    const int maxInitLevel = grid->getFinestInitializedLevel();

    std::vector<std::vector<SPtr<Block3D>>> blockVector(maxInitLevel + 1);
    for (int level = minInitLevel; level <= maxInitLevel; ++level) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);
    }

    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8>  cells;
    std::vector<std::string>   datanames;
    datanames.emplace_back("Winding");
    datanames.emplace_back("Level");
    std::vector<std::vector<double>> data(datanames.size());

    for (int level = minInitLevel; level <= maxInitLevel; ++level) {
        for (const auto &block : blockVector[level]) {
            if (!block)
                continue;
            auto kernel = block->getKernel();
            if (!kernel)
                continue;
            auto bcArray = kernel->getBCSet()->getBCArray();
            if (!bcArray)
                continue;

            const double deltaCell = static_cast<double>(grid->getDeltaX(block));
            double       faceEps   = 0.0;
#if defined(VF_HAS_FAST_WINDING)
            if (hasWinding)
                faceEps = std::max(eval.faceEpsBase,
                                   ub_math::greater(deltaCell, 0.0) ? deltaCell * 1.0e-6 : 0.0);
#endif

            int minX1 = 0, minX2 = 0, minX3 = 0;
            int maxX1 = static_cast<int>(bcArray->getNX1());
            int maxX2 = static_cast<int>(bcArray->getNX2());
            int maxX3 = static_cast<int>(bcArray->getNX3());

            CbArray3D<int> nodeNumbers(maxX1, maxX2, maxX3, -1);
            int nr = static_cast<int>(nodes.size());

            maxX1 -= 2;
            maxX2 -= 2;
            maxX3 -= 2;

            for (int ix3 = minX3; ix3 <= maxX3; ++ix3) {
                for (int ix2 = minX2; ix2 <= maxX2; ++ix2) {
                    for (int ix1 = minX1; ix1 <= maxX1; ++ix1) {
                        nodeNumbers(ix1, ix2, ix3) = nr++;
                        GbVector3D world = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                        nodes.emplace_back(static_cast<float>(world[0]),
                                           static_cast<float>(world[1]),
                                           static_cast<float>(world[2]));

                        float winding = 0.0f;
#if defined(VF_HAS_FAST_WINDING)
                        if (hasWinding) {
                            float computed = 0.0f;
                            if (eval.computeClampedIfInside(world[0], world[1], world[2], faceEps, computed))
                                winding = computed;
                        }
#endif
                        data[0].push_back(static_cast<double>(winding));
                        data[1].push_back(static_cast<double>(level));
                    }
                }
            }

            // build cells (octs)
            int SWB = 0, SEB = 0, NEB = 0, NWB = 0, SWT = 0, SET = 0, NET = 0, NWT = 0;
            maxX1 -= 1;
            maxX2 -= 1;
            maxX3 -= 1;
            for (int ix3 = minX3; ix3 <= maxX3; ++ix3) {
                for (int ix2 = minX2; ix2 <= maxX2; ++ix2) {
                    for (int ix1 = minX1; ix1 <= maxX1; ++ix1) {
                        if ((SWB = nodeNumbers(ix1, ix2, ix3)) >= 0 && (SEB = nodeNumbers(ix1 + 1, ix2, ix3)) >= 0 &&
                            (NEB = nodeNumbers(ix1 + 1, ix2 + 1, ix3)) >= 0 &&
                            (NWB = nodeNumbers(ix1, ix2 + 1, ix3)) >= 0 &&
                            (SWT = nodeNumbers(ix1, ix2, ix3 + 1)) >= 0 &&
                            (SET = nodeNumbers(ix1 + 1, ix2, ix3 + 1)) >= 0 &&
                            (NET = nodeNumbers(ix1 + 1, ix2 + 1, ix3 + 1)) >= 0 &&
                            (NWT = nodeNumbers(ix1, ix2 + 1, ix3 + 1)) >= 0) {
                            cells.emplace_back(static_cast<unsigned int>(SWB), static_cast<unsigned int>(SEB),
                                               static_cast<unsigned int>(NEB), static_cast<unsigned int>(NWB),
                                               static_cast<unsigned int>(SWT), static_cast<unsigned int>(SET),
                                               static_cast<unsigned int>(NET), static_cast<unsigned int>(NWT));
                        }
                    }
                }
            }
        }
    }

    const int istep         = static_cast<int>(step);
    std::string subfolder   = "winding" + ub_system::toString(istep);
    std::string pfilePath   = path_ + "/bc/" + subfolder;
    std::string cfilePath   = path_ + "/bc/winding_collection";
    std::string partPath    = pfilePath + "/winding" + ub_system::toString(gridRank) + "_" + ub_system::toString(istep);

    std::string partName = WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(
        partPath, nodes, cells, datanames, data);

    const std::size_t found = partName.find_last_of("/");
    std::string piece       = partName.substr(found + 1);
    piece                   = subfolder + "/" + piece;

    std::vector<std::string> cellDataNames;
    std::vector<std::string> pieces = comm ? comm->gather(piece) : std::vector<std::string>{ piece };

    if (!comm || comm->isRoot()) {
        std::string pname = WbWriterVtkXmlBinary::getInstance()->writeParallelFile(
            pfilePath, pieces, datanames, cellDataNames);

        const std::size_t foundName = pname.find_last_of("/");
        std::string       collectionPiece = pname.substr(foundName + 1);
        std::vector<std::string> filenames{ collectionPiece };

        if (ub_math::equal(static_cast<double>(step), scheduler->getMinBegin())) {
            WbWriterVtkXmlBinary::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else {
            WbWriterVtkXmlBinary::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }

        VF_LOG_INFO("GridWindingDiagnosticsSimulationObserver - winding volume written for step {}", istep);
    }
}

//! \}
