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
//! \addtogroup cpu_Interactors Interactors
//! \ingroup cpu_core core
//! \{
//! \author Hussein Alihussein
//! \note Helper interactor for CPU grid-winding diagnostics, introduced with assistance from GPT-5.1.
//=======================================================================================

#include "D3Q27GridWindingInteractor.h"

#include "BC.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "LBMKernel.h"
#include "Grid3D.h"

#include <basics/Timer/Timer.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/geometry3d/GbVector3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>

// #include <basics/geometry3d/winding/GeneralizedWindingNumber.h>
#include <basics/geometry3d/winding/GridWindingClassification.h>
#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>
#include <basics/geometry3d/winding/GridWindingWriting.h>

#include <logger/Logger.h>
#include <parallel/Communicator.h>
#if defined(VF_MPI)
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

void D3Q27GridWindingInteractor::initInteractor(const real &timeStep)
{
    auto gridPtr = grid.lock();
    if (!gridPtr)
        return;

    setQs(timeStep);
}

//////////////////////////////////////////////////////////////////////////

void D3Q27GridWindingInteractor::updateInteractor(const real &timestep)
{
    D3Q27Interactor::updateInteractor(timestep);
}

//////////////////////////////////////////////////////////////////////////

bool D3Q27GridWindingInteractor::setDifferencesToGbObject3D(const SPtr<Block3D> block)
{
    if (!block || block->isNotActive())
        return false;

    auto kernel = block->getKernel();
    if (!kernel)
        return false;

    auto bcSet = kernel->getBCSet();
    if (!bcSet)
        return false;

    auto bcArray = bcSet->getBCArray();
    if (!bcArray)
        return false;

    auto &bcNodes    = bcNodeIndicesMap[block];
    auto &solidNodes = solidNodeIndicesMap[block];
    bcNodes.clear();
    solidNodes.clear();

    bool oneEntryGotBC = false;

    const int nx1 = static_cast<int>(bcArray->getNX1());
    const int nx2 = static_cast<int>(bcArray->getNX2());
    const int nx3 = static_cast<int>(bcArray->getNX3());

    for (int iz = 0; iz < nx3; ++iz) {
        for (int iy = 0; iy < nx2; ++iy) {
            for (int ix = 0; ix < nx1; ++ix) {
                if (bcArray->isSolid(ix, iy, iz))
                    solidNodes.insert(UbTupleInt3(ix, iy, iz));
                if (bcArray->hasBC(ix, iy, iz)) {
                    bcNodes.insert({ ix, iy, iz });
                    oneEntryGotBC = true;
                }
            }
        }
    }

    return oneEntryGotBC;
}

//////////////////////////////////////////////////////////////////////////

void D3Q27GridWindingInteractor::setQs(const real &timeStep)
{
#if !defined(VF_HAS_FAST_WINDING)
    UB_THROW(UbException(UB_EXARGS, "D3Q27GridWindingInteractor requires VF_HAS_FAST_WINDING"));
#endif

    auto gridPtr = grid.lock();
    if (!gridPtr)
        return;

    auto surface = std::dynamic_pointer_cast<GbTriFaceMesh3D>(geoObject3D);
    if (!surface) {
        D3Q27Interactor::initInteractor(timeStep);
        return;
    }

    auto bcSet = findBCSet(gridPtr);
    if (!bcSet)
        return;

    vf::grid_winding::cpu::BoundaryProcessingConfig config;

    auto comm = vf::parallel::Communicator::getInstance();
    computeBoundaryData(surface, gridPtr, bcSet, BCs, static_cast<double>(timeStep), comm, std::string(), config);

    bcNodeIndicesMap.clear();
    solidNodeIndicesMap.clear();
    if (bcBlocks.empty())
        rebuildNodeIndexMaps();
    else
        updateBlocks();

    bool needTimeDependence = false;
    for (const auto &bc : BCs) {
        if (bc && bc->isTimeDependent())
            needTimeDependence = true;
    }
    if (needTimeDependence)
        this->setTimeDependent();
    else
        this->unsetTimeDependent();
}

//////////////////////////////////////////////////////////////////////////

void D3Q27GridWindingInteractor::computeBoundaryData(
    const SPtr<GbTriFaceMesh3D>            &surface,
    const SPtr<Grid3D>                     &grid,
    const SPtr<BCSet>                      &bcSet,
    const std::vector<SPtr<BC>>            &boundaryAdapters,
    double                                  timeStep,
    const SPtr<vf::parallel::Communicator> &comm,
    const std::string                      &basePath,
    const vf::grid_winding::cpu::BoundaryProcessingConfig &config)
{
    if (!surface || !grid || !bcSet)
        return;

    surface->calculateValues();

    auto effectiveConfig = config;
    // Let callers decide whether files are written; this helper
    // simply forwards the configuration.
    boundaryResult_ = vf::grid_winding::cpu::processBoundaryData(
        surface, grid, bcSet, boundaryAdapters, this, timeStep, comm, basePath, effectiveConfig);
}

//////////////////////////////////////////////////////////////////////////

void D3Q27GridWindingInteractor::rebuildNodeIndexMaps()
{
    bcNodeIndicesMap.clear();
    solidNodeIndicesMap.clear();

    auto gridPtr = grid.lock();
    if (!gridPtr)
        return;

    const int coarsestLevel = gridPtr->getCoarsestInitializedLevel();
    const int finestLevel   = gridPtr->getFinestInitializedLevel();

    for (int level = coarsestLevel; level <= finestLevel; ++level) {
        std::vector<SPtr<Block3D>> blocks;
        gridPtr->getBlocks(level, blocks);

        for (const auto &block : blocks)
            setDifferencesToGbObject3D(block);
    }
}

//////////////////////////////////////////////////////////////////////////

SPtr<BCSet> D3Q27GridWindingInteractor::findBCSet(const SPtr<Grid3D> &grid) const
{
    if (!grid)
        return {};

    const int coarsestLevel = grid->getCoarsestInitializedLevel();
    const int finestLevel   = grid->getFinestInitializedLevel();

    for (int level = coarsestLevel; level <= finestLevel; ++level) {
        std::vector<SPtr<Block3D>> blocks;
        grid->getBlocks(level, blocks);

        for (const auto &block : blocks) {
            if (!block)
                continue;
            auto kernel = block->getKernel();
            if (!kernel)
                continue;
            auto bcSet = kernel->getBCSet();
            if (bcSet)
                return bcSet;
        }
    }
    return {};
}

//////////////////////////////////////////////////////////////////////////

namespace vf::grid_winding::cpu
{
static QLineCollection collectQLines(const SPtr<Grid3D> &grid,
                                     const SPtr<vf::parallel::Communicator> &comm)
{
    QLineCollection collection;
    if (!grid)
        return collection;

    std::vector<double> coordsSend;
    std::vector<float>  qSend;
    std::vector<float>  qRawSend;
    std::vector<float>  deltaSend;
    std::vector<int>    dirSend;

    const int coarsest = grid->getCoarsestInitializedLevel();
    const int finest   = grid->getFinestInitializedLevel();

    for (int level = coarsest; level <= finest; ++level) {
        std::vector<SPtr<Block3D>> blocks;
        grid->getBlocks(level, blocks);
        for (const auto &block : blocks) {
            if (!block)
                continue;

            auto kernel = block->getKernel();
            if (!kernel)
                continue;

            auto bcSet = kernel->getBCSet();
            if (!bcSet)
                continue;

            auto bcArray = bcSet->getBCArray();
            if (!bcArray)
                continue;

            const int nx = static_cast<int>(bcArray->getNX1());
            const int ny = static_cast<int>(bcArray->getNX2());
            const int nz = static_cast<int>(bcArray->getNX3());

            const double dx = static_cast<double>(grid->getDeltaX(block));

            for (int iz = 0; iz < nz; ++iz) {
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        if (!bcArray->hasBC(ix, iy, iz))
                            continue;

                        auto bc = bcArray->getBC(ix, iy, iz);
                        if (!bc)
                            continue;

                        GbVector3D originVec = grid->getNodeCoordinates(block, ix, iy, iz);
                        const double ox      = originVec[0];
                        const double oy      = originVec[1];
                        const double oz      = originVec[2];

                        for (int dir = d3q27_system::FSTARTDIR; dir <= d3q27_system::FENDDIR; ++dir) {
                            const float rawQ = bc->getQ(dir);
                            if (ub_math::less(rawQ, 0.0f))
                                continue;

                            const int cx = d3q27_system::DX1[dir];
                            const int cy = d3q27_system::DX2[dir];
                            const int cz = d3q27_system::DX3[dir];
                            if (cx == 0 && cy == 0 && cz == 0)
                                continue;

                            const float clampedQ = std::clamp(rawQ, 0.0f, 1.0f);

                            const double hitX = ox + static_cast<double>(clampedQ) * dx * static_cast<double>(cx);
                            const double hitY = oy + static_cast<double>(clampedQ) * dx * static_cast<double>(cy);
                            const double hitZ = oz + static_cast<double>(clampedQ) * dx * static_cast<double>(cz);

                            coordsSend.push_back(ox);
                            coordsSend.push_back(oy);
                            coordsSend.push_back(oz);
                            coordsSend.push_back(hitX);
                            coordsSend.push_back(hitY);
                            coordsSend.push_back(hitZ);

                            dirSend.push_back(dir);
                            qSend.push_back(clampedQ);
                            qRawSend.push_back(rawQ);
                            deltaSend.push_back(static_cast<float>(dx));
                        }
                    }
                }
            }
        }
    }

    std::vector<double> coordsRecv;
    std::vector<int>    dirRecv;
    std::vector<float>  qRecv;
    std::vector<float>  qRawRecv;
    std::vector<float>  deltaRecv;
    if (comm) {
        comm->allGather(coordsSend, coordsRecv);
        comm->allGather(dirSend, dirRecv);
        comm->allGather(qSend, qRecv);
        comm->allGather(qRawSend, qRawRecv);
        comm->allGather(deltaSend, deltaRecv);
    } else {
        coordsRecv = coordsSend;
        dirRecv    = dirSend;
        qRecv      = qSend;
        qRawRecv   = qRawSend;
        deltaRecv  = deltaSend;
    }

    if (comm && !comm->isRoot())
        return collection;

    if (coordsRecv.empty() || dirRecv.empty() || qRecv.empty())
        return collection;

    const std::size_t lineCount = qRecv.size();
    if (coordsRecv.size() < lineCount * 6)
        return collection;

    collection = buildQLineCollection(coordsRecv, dirRecv, qRecv);

    return collection;
}

//////////////////////////////////////////////////////////////////////////

// Adapter that exposes Grid3D nodes to the shared fast-winding solid marker logic.
class CpuWindingSolidAccessor
{
public:
    using Vec3 = vf::grid_winding::Vec3;

    struct Node
    {
        Vec3 position;
        BCArray3D *bcArray{ nullptr };
        int ix{ 0 };
        int iy{ 0 };
        int iz{ 0 };
    };

    CpuWindingSolidAccessor(const SPtr<Grid3D> &grid, std::array<double, 3> meshMin, std::array<double, 3> meshMax)
        : grid_(grid), meshMin_(meshMin), meshMax_(meshMax)
    {
        if (grid_) {
            const int coarsestLevel = grid_->getCoarsestInitializedLevel();
            std::vector<SPtr<Block3D>> blocks;
            grid_->getBlocks(coarsestLevel, blocks);
            if (!blocks.empty())
                cellSizeHint_ = static_cast<double>(grid_->getDeltaX(blocks.front()));
        }
    }

    template <typename Fn>
    void forEachNode(Fn &&fn)
    {
        const double faceEpsBase = vf::grid_winding::computeFaceEpsilon(meshMin_, meshMax_, cellSizeHint_);

        if (!grid_)
            return;

        const int coarsestLevel = grid_->getCoarsestInitializedLevel();
        const int finestLevel   = grid_->getFinestInitializedLevel();

        for (int level = coarsestLevel; level <= finestLevel; ++level) {
            std::vector<SPtr<Block3D>> blocks;
            grid_->getBlocks(level, blocks);
            for (const auto &block : blocks) {
                auto kernel = block->getKernel();
                if (!kernel)
                    continue;

                auto bcSet = kernel->getBCSet();
                if (!bcSet)
                    continue;

                auto bcArray = bcSet->getBCArray();
                if (!bcArray)
                    continue;

                real minX1, minX2, minX3, maxX1, maxX2, maxX3;

                UbTupleDouble3 blockLengths = grid_->getBlockLengths(block);
                UbTupleDouble3 org          = grid_->getBlockWorldCoordinates(block);
                UbTupleDouble3 nodeOffset   = grid_->getNodeOffset(block);

                minX1 = val<1>(org) + val<1>(nodeOffset);
                minX2 = val<2>(org) + val<2>(nodeOffset);
                minX3 = val<3>(org) + val<3>(nodeOffset);
                maxX1 = val<1>(org) + val<1>(blockLengths) - val<1>(nodeOffset);
                maxX2 = val<2>(org) + val<2>(blockLengths) - val<2>(nodeOffset);
                maxX3 = val<3>(org) + val<3>(blockLengths) - val<3>(nodeOffset);

                const double deltaCell = static_cast<double>(grid_->getDeltaX(block));
                const double faceEps =
                    std::max(faceEpsBase, ub_math::greater(deltaCell, 0.0) ? deltaCell * 1.0e-6 : 0.0);

                if (ub_math::less(maxX1, meshMin_[0] - faceEps) ||
                    ub_math::greater(minX1, meshMax_[0] + faceEps) ||
                    ub_math::less(maxX2, meshMin_[1] - faceEps) ||
                    ub_math::greater(minX2, meshMax_[1] + faceEps) ||
                    ub_math::less(maxX3, meshMin_[2] - faceEps) ||
                    ub_math::greater(minX3, meshMax_[2] + faceEps)) {
                    continue;
                }

                const int NX = static_cast<int>(bcArray->getNX1()) - 1;
                const int NY = static_cast<int>(bcArray->getNX2()) - 1;
                const int NZ = static_cast<int>(bcArray->getNX3()) - 1;

                for (int k = 0; k <= NZ; ++k)
                    for (int j = 0; j <= NY; ++j)
                        for (int i = 0; i <= NX; ++i) {
                            GbVector3D coord = grid_->getNodeCoordinates(block, i, j, k);
                            Node node;
                            node.position = Vec3{ coord[0], coord[1], coord[2] };
                            node.bcArray  = bcArray.get();
                            node.ix       = i;
                            node.iy       = j;
                            node.iz       = k;
                            fn(node);
                        }
            }
        }
    }

    void setSolid(const Node &node, float /*windingValue*/) const
    {
        if (!node.bcArray)
            return;
        if (!node.bcArray->isSolid(node.ix, node.iy, node.iz)) {
            node.bcArray->setSolid(static_cast<std::size_t>(node.ix), static_cast<std::size_t>(node.iy),
                                   static_cast<std::size_t>(node.iz));
        }
    }

    [[nodiscard]] double cellSizeHint() const
    {
        return cellSizeHint_;
    }

private:
    SPtr<Grid3D> grid_;
    std::array<double, 3> meshMin_;
    std::array<double, 3> meshMax_;
    double cellSizeHint_{ 0.0 };
};

// void markSolidsWithWinding(const SPtr<Grid3D> &grid, const vf::geometry::GeneralizedWindingNumber &winding,
//                            const std::array<double, 3> &meshMin, const std::array<double, 3> &meshMax,
//                            double threshold, double tolerance)
// {
//     CpuWindingSolidAccessor accessor(grid, meshMin, meshMax);
//     const double interiorBias = 2.0 * tolerance;
//     vf::grid_winding::markSolidsWithWinding(accessor, winding, meshMin, meshMax, threshold + interiorBias, 0.0);
// }

#if defined(VF_HAS_FAST_WINDING)
void markSolidsWithFastWinding(const SPtr<GbTriFaceMesh3D> &surface, const SPtr<Grid3D> &grid,
                               float accuracyScale, float threshold, float tolerance)
{
    if (!surface || !grid)
        return;

    const std::array<double, 3> meshMin{ surface->getX1Minimum(), surface->getX2Minimum(), surface->getX3Minimum() };
    const std::array<double, 3> meshMax{ surface->getX1Maximum(), surface->getX2Maximum(), surface->getX3Maximum() };

    CpuWindingSolidAccessor accessor(grid, meshMin, meshMax);
    vf::grid_winding::markSolidsWithFastWinding(accessor, *surface, accuracyScale, threshold, tolerance);
}

//////////////////////////////////////////////////////////////////////////

#endif

} // namespace vf::grid_winding::cpu

//////////////////////////////////////////////////////////////////////////

namespace vf::grid_winding::cpu
{
using vf::grid_winding::SubgridDistanceStats;
using vf::grid_winding::Vec3;
using vf::grid_winding::processBoundaryLink;

namespace detail
{
inline void ensureBoundaryConditions(const SPtr<BCArray3D> &bcArray, int ix, int iy, int iz,
                                     SPtr<BoundaryConditions> &outBc)
{
    outBc = bcArray->getBC(ix, iy, iz);
    if (!outBc) {
        outBc = std::make_shared<BoundaryConditions>();
        bcArray->setBC(ix, iy, iz, outBc);
    }
}

struct BoundaryLink
{
    int ix{ 0 };
    int iy{ 0 };
    int iz{ 0 };
    int dir{ 0 };
    bool fluidLocal{ false };
    bool solidLocal{ false };
};

struct RemoteBoundaryLinkMessage
{
    int blockX1{ 0 };
    int blockX2{ 0 };
    int blockX3{ 0 };
    int level{ 0 };
    int dir{ 0 };
    int ix{ 0 };
    int iy{ 0 };
    int iz{ 0 };
};

static_assert(std::is_trivially_copyable_v<RemoteBoundaryLinkMessage>,
              "RemoteBoundaryLinkMessage must be trivially copyable for MPI exchange");

} // namespace detail

inline void computeSubgridDistancesStandalone(const SPtr<GbTriFaceMesh3D> &surface, const SPtr<Grid3D> &grid,
                                              const SPtr<BCSet> &bcSet,
                                              const std::vector<SPtr<BC>> &boundaryAdapters,
                                              const ::D3Q27Interactor *adapterInteractor,
                                              double timeStep,
                                              SubgridDistanceStats *stats)
{
    using namespace detail;

    if (!surface || !grid || !bcSet)
        return;

    auto *triangles = surface->getTriangles();
    auto *nodes = surface->getNodes();
    if (!nodes || !triangles || triangles->empty())
        return;

    surface->calculateValues();

    const std::array<double, 3> meshMin{ surface->getX1Minimum(), surface->getX2Minimum(), surface->getX3Minimum() };
    const std::array<double, 3> meshMax{ surface->getX1Maximum(), surface->getX2Maximum(), surface->getX3Maximum() };

    std::vector<GbTriFaceMesh3D *> meshes{ surface.get() };
    auto kdTrees = vf::grid_winding::ensureKdTrees(meshes);

    const bool hasAdapters = (adapterInteractor != nullptr && !boundaryAdapters.empty());
    if (hasAdapters) {
        for (const auto &adapter : boundaryAdapters) {
            if (adapter)
                adapter->init(adapterInteractor, timeStep);
        }
    }

    auto communicator = vf::parallel::Communicator::getInstance();
    int currentRank   = communicator ? communicator->getProcessID() : 0;
    int worldSize     = communicator ? communicator->getNumberOfProcesses() : 1;
    if (worldSize <= 0)
        worldSize = 1;

    const bool isRootRank = !communicator || communicator->isRoot();
    if (isRootRank && meshes.size() > 0 && kdTrees.empty()) {
        VF_LOG_WARNING("No kd-tree available for grid-winding; falling back to triangle list intersections.");
    }

#if defined(VF_MPI)
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    if (communicator) {
        if (void *native = communicator->getNativeCommunicator()) {
            mpiComm = *static_cast<MPI_Comm *>(native);
        }
        currentRank = communicator->getProcessID();
        worldSize   = communicator->getNumberOfProcesses();
        if (worldSize <= 0)
            worldSize = 1;
    } else {
        MPI_Comm_rank(mpiComm, &currentRank);
        MPI_Comm_size(mpiComm, &worldSize);
        if (worldSize <= 0)
            worldSize = 1;
    }
#endif

    std::vector<std::vector<detail::RemoteBoundaryLinkMessage>> remoteMessages(
        static_cast<std::size_t>(worldSize));

    struct PromotedSolid
    {
        SPtr<Block3D> block;
        int           ix{ 0 };
        int           iy{ 0 };
        int           iz{ 0 };
    };

    std::vector<PromotedSolid> promotedSolids;
    promotedSolids.reserve(128);

    struct CpuKdAdapter
    {
        BCArray3D *bcArray{ nullptr };
        const std::vector<SPtr<BC>> *boundaryAdapters{ nullptr };
        const ::D3Q27Interactor     *adapterInteractor{ nullptr };
        bool                         hasAdapters{ false };
        SubgridDistanceStats        *stats{ nullptr };
        SPtr<Block3D>                ownerBlock;
        std::vector<PromotedSolid>  *promotedSolids{ nullptr };
        int                          ix{ 0 };
        int                          iy{ 0 };
        int                          iz{ 0 };
        int                          dir{ 0 };
        bool                         solidLocal{ false };
        Vec3                         origin{};

        bool promoteSolid(const vf::grid_winding::BoundaryLink &) const
        {
            if (!bcArray)
                return false;
            if (!bcArray->isFluid(static_cast<std::size_t>(ix),
                                  static_cast<std::size_t>(iy),
                                  static_cast<std::size_t>(iz)))
                return false;
            bcArray->setSolid(static_cast<std::size_t>(ix),
                              static_cast<std::size_t>(iy),
                              static_cast<std::size_t>(iz));
            return true;
        }

        double existingQ(const vf::grid_winding::BoundaryLink &) const
        {
            if (!bcArray)
                return -1.0;
            auto bc = bcArray->getBC(ix, iy, iz);
            if (!bc)
                return -1.0;
            return static_cast<double>(bc->getQ(dir));
        }

        void setQ(const vf::grid_winding::BoundaryLink &, double q, int) const
        {
            if (!bcArray)
                return;
            if (bcArray->isSolid(static_cast<std::size_t>(ix),
                                 static_cast<std::size_t>(iy),
                                 static_cast<std::size_t>(iz)))
                return;

            // if (ub_math::zero(q)) {
            //     bcArray->setSolid(static_cast<std::size_t>(ix),
            //                       static_cast<std::size_t>(iy),
            //                       static_cast<std::size_t>(iz));
            //     if (stats)
            //         ++stats->promotedZeroQNodes;
            //     if (promotedSolids && ownerBlock)
            //         promotedSolids->push_back(PromotedSolid{ ownerBlock, ix, iy, iz });
            //     return;
            // }

            auto bc = bcArray->getBC(ix, iy, iz);
            if (!bc) {
                bc = std::make_shared<BoundaryConditions>();
                bcArray->setBC(ix, iy, iz, bc);
            }
            bc->setQ(static_cast<float>(q), dir);

            if (hasAdapters && boundaryAdapters && adapterInteractor) {
                const double ox = origin.X1();
                const double oy = origin.X2();
                const double oz = origin.X3();

                for (auto it = boundaryAdapters->rbegin(); it != boundaryAdapters->rend(); ++it) {
                    if (*it)
                        (*it)->adaptBCForDirection(*adapterInteractor, bc, ox, oy, oz, q, dir);
                }
                for (auto it = boundaryAdapters->rbegin(); it != boundaryAdapters->rend(); ++it) {
                    if (*it)
                        (*it)->adaptBC(*adapterInteractor, bc, ox, oy, oz);
                }
            }
        }

        void markSurface(const vf::grid_winding::BoundaryLink &, double, int) const {}

        std::pair<bool, double> windingNumber(const vf::grid_winding::BoundaryLink &) const
        {
            return { false, 0.0 };
        }

        bool handleMissing(const vf::grid_winding::BoundaryLink &,
                           SubgridDistanceStats::MissingLink &) const
        {
            if (!bcArray)
                return false;

            auto bc = bcArray->getBC(ix, iy, iz);
            if (!bc) {
                bc = std::make_shared<BoundaryConditions>();
                bcArray->setBC(ix, iy, iz, bc);
            }

            const double fallbackQ = vf::grid_winding::defaultMissingQ();
            bc->setQ(static_cast<float>(fallbackQ), dir);

            if (hasAdapters && boundaryAdapters && adapterInteractor) {
                const double ox = origin.X1();
                const double oy = origin.X2();
                const double oz = origin.X3();
                for (auto it = boundaryAdapters->rbegin(); it != boundaryAdapters->rend(); ++it) {
                    if (*it)
                        (*it)->adaptBCForDirection(*adapterInteractor, bc, ox, oy, oz,
                                                   fallbackQ, dir);
                }
                for (auto it = boundaryAdapters->rbegin(); it != boundaryAdapters->rend(); ++it) {
                    if (*it)
                        (*it)->adaptBC(*adapterInteractor, bc, ox, oy, oz);
                }
            }
            return true;
        }
    };

    CpuKdAdapter adapter;

    auto processLinkForBlock = [&](const SPtr<Block3D> &ownerBlock, int ix, int iy, int iz, int dir,
                                   SubgridDistanceStats *linkStats, bool allowGhost) {
        if (!ownerBlock)
            return;

        auto kernel = ownerBlock->getKernel();
        if (!kernel)
            return;

        auto bcSet = kernel->getBCSet();
        if (!bcSet)
            return;

        auto bcArray = bcSet->getBCArray();
        if (!bcArray)
            return;

        // if (!allowGhost && !kernel->isInsideOfDomain(ix, iy, iz)) {
        //     if (linkStats)
        //         ++linkStats->skippedNonLocalLinks;
        //     return;
        // }

        if (linkStats)
            ++linkStats->totalLinks;

        if (!bcArray->isFluid(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy),
                              static_cast<std::size_t>(iz))) {
            if (bcArray->isSolid(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy),
                                 static_cast<std::size_t>(iz)))
                return;
            bcArray->setFluid(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy),
                              static_cast<std::size_t>(iz));
        }

        const double dxLocal = static_cast<double>(grid->getDeltaX(ownerBlock));
        GbVector3D originVec = grid->getNodeCoordinates(ownerBlock, ix, iy, iz);
        Vec3       origin{ originVec[0], originVec[1], originVec[2] };

        const int cx = d3q27_system::DX1[dir];
        const int cy = d3q27_system::DX2[dir];
        const int cz = d3q27_system::DX3[dir];

        const int  solidX    = ix + cx;
        const int  solidY    = iy + cy;
        const int  solidZ    = iz + cz;
        const bool solidLocal = kernel->isInsideOfDomain(solidX, solidY, solidZ);

        adapter.bcArray          = bcArray.get();
        adapter.boundaryAdapters = &boundaryAdapters;
        adapter.adapterInteractor = adapterInteractor;
        adapter.hasAdapters      = hasAdapters;
        adapter.stats            = linkStats;
        adapter.ownerBlock       = ownerBlock;
        adapter.promotedSolids   = &promotedSolids;
        adapter.ix               = ix;
        adapter.iy               = iy;
        adapter.iz               = iz;
        adapter.dir              = dir;
        adapter.solidLocal       = solidLocal;
        adapter.origin           = origin;

        vf::grid_winding::BoundaryLink link;
        link.origin    = origin;
        link.neighbour = Vec3{ origin.X1() + dxLocal * static_cast<double>(cx),
                               origin.X2() + dxLocal * static_cast<double>(cy),
                               origin.X3() + dxLocal * static_cast<double>(cz) };
        link.direction = dir;

        const double faceEps = vf::grid_winding::computeFaceEpsilon(meshMin, meshMax, dxLocal);
        auto insideGeometry = [&](const Vec3 &p) {
            return ub_math::greaterEqual(p.X1(), meshMin[0] - faceEps) &&
                   ub_math::lessEqual(p.X1(), meshMax[0] + faceEps) &&
                   ub_math::greaterEqual(p.X2(), meshMin[1] - faceEps) &&
                   ub_math::lessEqual(p.X2(), meshMax[1] + faceEps) &&
                   ub_math::greaterEqual(p.X3(), meshMin[2] - faceEps) &&
                   ub_math::lessEqual(p.X3(), meshMax[2] + faceEps);
        };
        if (!insideGeometry(link.origin) && !insideGeometry(link.neighbour)) {
            return;
        }

        if (!kdTrees.empty()) {
            vf::grid_winding::QComputationConfig cfg;
            vf::grid_winding::processBoundaryLinkKdTree(kdTrees, adapter, link, linkStats, cfg);
        } else {
            auto promoteSolid = [&]() -> bool {
                return adapter.promoteSolid(link);
            };

            auto existingQ = [&]() -> double {
                return adapter.existingQ(link);
            };

            auto setQ = [&](double q, int patchIndex) {
                adapter.setQ(link, q, patchIndex);
            };

            auto markSurface = [&](double q, int patchIndex) {
                adapter.markSurface(link, q, patchIndex);
            };

            auto getWindingValue = [&]() -> std::pair<bool, double> {
                return adapter.windingNumber(link);
            };

            auto handleMissing = [&](int /*missingDir*/, const Vec3 &, const Vec3 &,
                                     SubgridDistanceStats::MissingLink &missing) -> bool {
                return adapter.handleMissing(link, missing);
            };

            processBoundaryLink(meshes, origin, dxLocal, cx, cy, cz, dir, promoteSolid, existingQ,
                                setQ, markSurface, getWindingValue, handleMissing, linkStats);
        }
    };

    const int coarsestLevel = grid->getCoarsestInitializedLevel();
    const int finestLevel   = grid->getFinestInitializedLevel();

    for (int level = coarsestLevel; level <= finestLevel; ++level) {
        std::vector<SPtr<Block3D>> blocks;
        grid->getBlocks(level, blocks);
        for (const auto &block : blocks) {
            auto kernel = block->getKernel();
            if (!kernel)
                continue;

            auto blockBcSet = kernel->getBCSet();
            if (!blockBcSet)
                continue;
            auto bcArray = blockBcSet->getBCArray();
            if (!bcArray)
                continue;

            const int nx = static_cast<int>(bcArray->getNX1());
            const int ny = static_cast<int>(bcArray->getNX2());
            const int nz = static_cast<int>(bcArray->getNX3());

            std::vector<BoundaryLink> boundaryLinks;
            boundaryLinks.reserve(static_cast<std::size_t>(nx * ny));

            // Collect boundary links by scanning solid cells for adjacent fluid neighbours.
            for (int iz = 0; iz < nz; ++iz) {
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        if (!bcArray->isSolid(ix, iy, iz))
                            continue;

                        // const bool solidLocal = kernel->isInsideOfDomain(ix, iy, iz);
                        bool solidLocal = true;

                        for (int dir = d3q27_system::FSTARTDIR; dir <= d3q27_system::FENDDIR; ++dir) {
                            const int cx = d3q27_system::DX1[dir];
                            const int cy = d3q27_system::DX2[dir];
                            const int cz = d3q27_system::DX3[dir];
                            if (cx == 0 && cy == 0 && cz == 0)
                                continue;

                            const int fluidX = ix + cx;
                            const int fluidY = iy + cy;
                            const int fluidZ = iz + cz;
                            if (fluidX < 0 || fluidY < 0 || fluidZ < 0 || fluidX >= nx || fluidY >= ny ||
                                fluidZ >= nz)
                                continue;

                            if (bcArray->isSolid(static_cast<std::size_t>(fluidX), static_cast<std::size_t>(fluidY),
                                                 static_cast<std::size_t>(fluidZ)))
                                continue;

                            // const bool fluidLocal = kernel->isInsideOfDomain(fluidX, fluidY, fluidZ);
                            bool fluidLocal = true;
                            if (!solidLocal && !fluidLocal)
                                continue;

                            const int fluidDir = d3q27_system::INVDIR[dir];
                            if (fluidDir < d3q27_system::FSTARTDIR || fluidDir > d3q27_system::FENDDIR)
                                continue;

                            boundaryLinks.push_back(
                                BoundaryLink{ fluidX, fluidY, fluidZ, fluidDir, fluidLocal, solidLocal });
                        }
                    }
                }
            }

            for (int iz = 0; iz < nz; ++iz) {
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        // if (!kernel->isInsideOfDomain(ix, iy, iz))
                        //     continue;
                        if (bcArray->isSolid(ix, iy, iz))
                            continue;

                        const bool originFluid = bcArray->isFluid(ix, iy, iz);
                        const bool originUndefined = bcArray->isUndefined(ix, iy, iz);
                        if (!originFluid && !originUndefined)
                            continue;

                        for (int dir = d3q27_system::FSTARTDIR; dir <= d3q27_system::FENDDIR; ++dir) {
                            const int cx = d3q27_system::DX1[dir];
                            const int cy = d3q27_system::DX2[dir];
                            const int cz = d3q27_system::DX3[dir];
                            if (cx == 0 && cy == 0 && cz == 0)
                                continue;

                            const int solidX = ix + cx;
                            const int solidY = iy + cy;
                            const int solidZ = iz + cz;
                            if (solidX < 0 || solidY < 0 || solidZ < 0 || solidX >= nx || solidY >= ny ||
                                solidZ >= nz)
                                continue;

                            if (!bcArray->isSolid(static_cast<std::size_t>(solidX), static_cast<std::size_t>(solidY),
                                                  static_cast<std::size_t>(solidZ)))
                                continue;

                            // const bool solidLocal = kernel->isInsideOfDomain(solidX, solidY, solidZ);
                            const bool solidLocal = true;

                            boundaryLinks.push_back(
                                BoundaryLink{ ix, iy, iz, dir, true, solidLocal });
                        }
                    }
                }
            }

            if (boundaryLinks.empty())
                continue;

            std::sort(boundaryLinks.begin(), boundaryLinks.end(),
                      [](const BoundaryLink &lhs, const BoundaryLink &rhs) {
                          if (lhs.ix != rhs.ix)
                              return lhs.ix < rhs.ix;
                          if (lhs.iy != rhs.iy)
                              return lhs.iy < rhs.iy;
                          if (lhs.iz != rhs.iz)
                              return lhs.iz < rhs.iz;
                          return lhs.dir < rhs.dir;
                      });

            std::vector<BoundaryLink> uniqueLinks;
            uniqueLinks.reserve(boundaryLinks.size());
            for (const auto &link : boundaryLinks) {
                if (!uniqueLinks.empty() && uniqueLinks.back().ix == link.ix && uniqueLinks.back().iy == link.iy &&
                    uniqueLinks.back().iz == link.iz && uniqueLinks.back().dir == link.dir) {
                    auto &entry = uniqueLinks.back();
                    entry.fluidLocal = entry.fluidLocal || link.fluidLocal;
                    entry.solidLocal = entry.solidLocal || link.solidLocal;
                } else {
                    uniqueLinks.push_back(link);
                }
            }
            boundaryLinks.swap(uniqueLinks);

            for (const auto &link : boundaryLinks) {
                if (link.fluidLocal) {
                    processLinkForBlock(block, link.ix, link.iy, link.iz, link.dir, stats, false);
                    continue;
                }

                GbVector3D originVec = grid->getNodeCoordinates(block, link.ix, link.iy, link.iz);
                UbTupleInt3 ownerIndices =
                    grid->getBlockIndexes(originVec[0], originVec[1], originVec[2], block->getLevel());
                SPtr<Block3D> ownerBlock =
                    grid->getBlock(val<1>(ownerIndices), val<2>(ownerIndices), val<3>(ownerIndices), block->getLevel());
                if (!ownerBlock)
                    continue;

                UbTupleInt3 nodeIdx = grid->getNodeIndexes(ownerBlock, originVec[0], originVec[1], originVec[2]);
                const int ownerRank = ownerBlock->getRank();

                if (ownerRank == currentRank) {
                    processLinkForBlock(ownerBlock, val<1>(nodeIdx), val<2>(nodeIdx), val<3>(nodeIdx), link.dir,
                                        stats, false);
                } else if (ownerRank >= 0 && ownerRank < worldSize) {
                    if (stats)
                        ++stats->skippedNonLocalLinks;
                    detail::RemoteBoundaryLinkMessage msg;
                    msg.blockX1 = ownerBlock->getX1();
                    msg.blockX2 = ownerBlock->getX2();
                    msg.blockX3 = ownerBlock->getX3();
                    msg.level   = ownerBlock->getLevel();
                    msg.dir     = link.dir;
                    msg.ix      = val<1>(nodeIdx);
                    msg.iy      = val<2>(nodeIdx);
                    msg.iz      = val<3>(nodeIdx);
                    remoteMessages[static_cast<std::size_t>(ownerRank)].push_back(msg);
                }
            }
        }

    }

#if defined(VF_MPI)
    if (worldSize > 1) {
        std::vector<int> sendCounts(static_cast<std::size_t>(worldSize), 0);
        std::vector<int> sendDispls(static_cast<std::size_t>(worldSize), 0);
        std::vector<detail::RemoteBoundaryLinkMessage> sendBuffer;
        int offset = 0;
        for (int dest = 0; dest < worldSize; ++dest) {
            if (dest == currentRank) {
                sendDispls[static_cast<std::size_t>(dest)] = offset;
                continue;
            }
            auto &entries = remoteMessages[static_cast<std::size_t>(dest)];
            sendCounts[static_cast<std::size_t>(dest)] = static_cast<int>(entries.size());
            sendDispls[static_cast<std::size_t>(dest)] = offset;
            offset += sendCounts[static_cast<std::size_t>(dest)];
            sendBuffer.insert(sendBuffer.end(), entries.begin(), entries.end());
        }

        std::vector<int> recvCounts(static_cast<std::size_t>(worldSize), 0);
        MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, mpiComm);

        std::vector<int> recvDispls(static_cast<std::size_t>(worldSize), 0);
        int totalRecv = 0;
        for (int i = 0; i < worldSize; ++i) {
            recvDispls[static_cast<std::size_t>(i)] = totalRecv;
            totalRecv += recvCounts[static_cast<std::size_t>(i)];
        }
        std::vector<detail::RemoteBoundaryLinkMessage> recvBuffer(static_cast<std::size_t>(totalRecv));

        const int messageBytes = static_cast<int>(sizeof(detail::RemoteBoundaryLinkMessage));
        std::vector<int> sendCountsBytes(static_cast<std::size_t>(worldSize), 0);
        std::vector<int> sendDisplsBytes(static_cast<std::size_t>(worldSize), 0);
        std::vector<int> recvCountsBytes(static_cast<std::size_t>(worldSize), 0);
        std::vector<int> recvDisplsBytes(static_cast<std::size_t>(worldSize), 0);

        for (int i = 0; i < worldSize; ++i) {
            sendCountsBytes[static_cast<std::size_t>(i)] = sendCounts[static_cast<std::size_t>(i)] * messageBytes;
            sendDisplsBytes[static_cast<std::size_t>(i)] = sendDispls[static_cast<std::size_t>(i)] * messageBytes;
            recvCountsBytes[static_cast<std::size_t>(i)] = recvCounts[static_cast<std::size_t>(i)] * messageBytes;
            recvDisplsBytes[static_cast<std::size_t>(i)] = recvDispls[static_cast<std::size_t>(i)] * messageBytes;
        }

        MPI_Alltoallv(sendBuffer.data(), sendCountsBytes.data(), sendDisplsBytes.data(), MPI_BYTE, recvBuffer.data(),
                      recvCountsBytes.data(), recvDisplsBytes.data(), MPI_BYTE, mpiComm);

        for (const auto &msg : recvBuffer) {
            SPtr<Block3D> ownerBlock = grid->getBlock(msg.blockX1, msg.blockX2, msg.blockX3, msg.level);
            if (!ownerBlock)
                continue;
            if (ownerBlock->getRank() != currentRank)
                continue;
            processLinkForBlock(ownerBlock, msg.ix, msg.iy, msg.iz, msg.dir, stats, false);
        }
    }
#else
    (void)remoteMessages;
    (void)currentRank;
#endif

    if (!promotedSolids.empty()) {
        const std::size_t promotedCount = promotedSolids.size();
        for (std::size_t idx = 0; idx < promotedCount; ++idx) {
            const auto &promoted = promotedSolids[idx];
            auto block = promoted.block;
            if (!block)
                continue;

            auto kernel = block->getKernel();
            if (!kernel)
                continue;

            auto blockBcSet = kernel->getBCSet();
            if (!blockBcSet)
                continue;
            auto bcArray = blockBcSet->getBCArray();
            if (!bcArray)
                continue;

            const int nx = static_cast<int>(bcArray->getNX1());
            const int ny = static_cast<int>(bcArray->getNX2());
            const int nz = static_cast<int>(bcArray->getNX3());

            for (int dir = d3q27_system::FSTARTDIR; dir <= d3q27_system::FENDDIR; ++dir) {
                const int cx = d3q27_system::DX1[dir];
                const int cy = d3q27_system::DX2[dir];
                const int cz = d3q27_system::DX3[dir];
                if (cx == 0 && cy == 0 && cz == 0)
                    continue;

                const int fluidX = promoted.ix + cx;
                const int fluidY = promoted.iy + cy;
                const int fluidZ = promoted.iz + cz;
                if (fluidX < 0 || fluidY < 0 || fluidZ < 0 || fluidX >= nx || fluidY >= ny || fluidZ >= nz)
                    continue;

                if (bcArray->isSolid(static_cast<std::size_t>(fluidX), static_cast<std::size_t>(fluidY),
                                     static_cast<std::size_t>(fluidZ)))
                    continue;

                const int fluidDir = d3q27_system::INVDIR[dir];
                if (fluidDir < d3q27_system::FSTARTDIR || fluidDir > d3q27_system::FENDDIR)
                    continue;

                processLinkForBlock(block, fluidX, fluidY, fluidZ, fluidDir, stats, false);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////

} // namespace vf::grid_winding::cpu

//////////////////////////////////////////////////////////////////////////

namespace vf::grid_winding::cpu
{
BoundaryProcessingResult processBoundaryData(
    const SPtr<GbTriFaceMesh3D>               &surface,
    const SPtr<Grid3D>                        &grid,
    const SPtr<BCSet>                         &bcSet,
    const std::vector<SPtr<BC>>               &boundaryAdapters,
    const ::D3Q27Interactor                   *adapterInteractor,
    double                                     timeStep,
    const SPtr<vf::parallel::Communicator>    &comm,
    const std::string                         &basePath,
    const BoundaryProcessingConfig            &config)
{
    BoundaryProcessingResult result;

    if (!surface || !grid || !bcSet)
        return result;

    int rank = 0;
    if (comm) {
        rank = comm->getProcessID();
    }

#if defined(VF_HAS_FAST_WINDING)
    {
        vf::basics::Timer windingTimer;
        windingTimer.start();
        // Use shared defaults for accuracy/tolerance and keep the threshold configurable.
        vf::grid_winding::cpu::markSolidsWithFastWinding(
            surface, grid,
            vf::grid_winding::FastWindingDefaultAccuracyScale,
            config.fastWindingThreshold,
            vf::grid_winding::FastWindingDefaultTolerance);
        VF_LOG_INFO("CPU fast-winding classification completed on rank {} in {} sec",
                    rank,
                    windingTimer.getCurrentRuntimeInSeconds());
    }
#endif

    vf::basics::Timer qTimer;
    qTimer.start();

    vf::grid_winding::cpu::computeSubgridDistancesStandalone(
        surface, grid, bcSet, boundaryAdapters, adapterInteractor, timeStep, &result.stats);

    VF_LOG_INFO("CPU subgrid link computation finished on rank {} in {} sec "
                "(hits={}, total={}, promoted={}, skipped={})",
                rank,
                qTimer.getCurrentRuntimeInSeconds(),
                result.stats.hitLinks,
                result.stats.totalLinks,
                result.stats.promotedFluidNodes,
                result.stats.skippedNonLocalLinks);

    if (comm) {
        std::vector<double> promotedTotals{ static_cast<double>(result.stats.promotedFluidNodes) };
        std::vector<double> promotedZeroQ{ static_cast<double>(result.stats.promotedZeroQNodes) };
        comm->allReduceSum(promotedTotals);
        comm->allReduceSum(promotedZeroQ);
        if (comm->isRoot()) {
            VF_LOG_INFO("CPU total promotions (solid) total: {}", static_cast<uint64_t>(promotedTotals.front()));
            VF_LOG_INFO("CPU zero-q promotions (solid) total: {}", static_cast<uint64_t>(promotedZeroQ.front()));
        }
    } else {
        VF_LOG_INFO("CPU total promotions (solid) total: {}", result.stats.promotedFluidNodes);
        VF_LOG_INFO("CPU zero-q promotions (solid) total: {}", result.stats.promotedZeroQNodes);
    }

    const auto &stats = result.stats;
    if (stats.totalLinks == stats.hitLinks) {
        VF_LOG_INFO("All CPU boundary links intersected the geometry on rank {}.", rank);
    } else {
        VF_LOG_WARNING("Some CPU boundary links did not intersect the geometry on rank {}.", rank);
        if (stats.skippedNonLocalLinks > 0) {
           // VF_LOG_INFO("Skipped {} ghost boundary links that were not owned by rank {}.",
            //             stats.skippedNonLocalLinks, rank);
        }
    }

    if (config.writeMissingLinks) {
        vf::basics::Timer missingTimer;
        missingTimer.start();
        result.missingLinks = vf::grid_winding::collectMissingLinks(result.stats, comm);

        if (comm && comm->isRoot() && !basePath.empty()) {
            VF_LOG_INFO("CPU missing-link aggregation completed in {} sec",
                        missingTimer.getCurrentRuntimeInSeconds());
        }
    } else if (comm && comm->isRoot()) {
        VF_LOG_INFO("CPU missing-link diagnostics disabled; skipping aggregation.");
    }

    if (config.writeQLines) {
        vf::basics::Timer qCollectTimer;
        qCollectTimer.start();
        result.qLines = vf::grid_winding::cpu::collectQLines(
            grid, comm);

        if (comm && comm->isRoot() && !basePath.empty()) {
            VF_LOG_INFO("CPU Q-line collection completed in {} sec",
                        qCollectTimer.getCurrentRuntimeInSeconds());
        }
    } else if (comm && comm->isRoot()) {
        VF_LOG_INFO("CPU Q-line diagnostics disabled; skipping collection.");
    }

    return result;
}
} // namespace vf::grid_winding::cpu

//! \}
