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
//! \addtogroup gpu_geometries geometries
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "TriangularMeshStrategy.h"

#include "Timer/Timer.h"

#include "basics/geometry3d/GbTriFaceMesh3D.h"

#include "geometries/Triangle/Triangle.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "grid/GridImp.h"
#include "grid/NodeValues.h"

#include "utilities/GridWindingBoundaryGpuCore.h"
#include "utilities/GridWindingSolidMarkerGpu.h"
#include <logger/Logger.h>

namespace vf::gpu {

void TriangularMeshDiscretizationStrategy::appendFastWindingQSurfaces(
    GridImp* /*grid*/, std::vector<SPtr<GbTriFaceMesh3D>>& /*surfaces*/) const
{
}

void TriangularMeshDiscretizationStrategy::computeFastWindingQs(
    GridImp* /*grid*/, const std::vector<SPtr<GbTriFaceMesh3D>>& /*surfaces*/) const
{
}


void PointInObjectDiscretizationStrategy::doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
{
    triangularMesh->generateGbTriFaceMesh3D();

    VF_LOG_INFO("Start Point-In-Object Test");

    // trigger the GbTriFaceMesh3D to generate a kd-tree
    triangularMesh->getGbTriFaceMesh3D()->isPointInGbObject3D(0.0, 0.0, 0.0);

    vf::basics::Timer timer;
    timer.start();

    real outputTime = 60.0;
    
#pragma omp parallel for
    for (int index = 0; index < (int)grid->getSize(); index++)
    {
        if( grid->getFieldEntry(index) == InnerType ) continue;

        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);

        if (triangularMesh->getGbTriFaceMesh3D()->isPointInGbObject3D(x, y, z))
            grid->setNodeTo(index, InnerType);
        //else
        //    grid->setNodeTo(i, OuterType);

        if( timer.getCurrentRuntimeInSeconds() > outputTime ){
            VF_LOG_INFO("    {} / {} nodes tested", index, grid->getSize());
            timer.start();
        }
    }
    VF_LOG_INFO("Done Point-In-Object Test");
}


void RayCastingDiscretizationStrategy::doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char InnerType, char OuterType)
{
        auto mesh = triangularMesh->getGbTriFaceMesh3D();

        const real minXExact = triangularMesh->minmax.minX;
        const real minYExact = triangularMesh->minmax.minY;
        const real minZExact = triangularMesh->minmax.minZ;

        const real maxXExact = triangularMesh->minmax.maxX;
        const real maxYExact = triangularMesh->minmax.maxY;
        const real maxZExact = triangularMesh->minmax.maxZ;

        const auto min = grid->getMinimumOnNode(Vertex(minXExact, minYExact, minZExact));

        const real minX = min.x;
        const real minY = min.y;
        const real minZ = min.z;

        const auto max = grid->getMaximumOnNode(Vertex(maxXExact, maxYExact, maxZExact));

        const real maxX = max.x;
        const real maxY = max.y;
        const real maxZ = max.z;


        real x, y, z;
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (x = minX; x <= maxX; x += grid->getDelta())
                {
                    grid->setNodeTo(grid->transCoordToIndex(x, y, z), InnerType);
                }
            }
        }

        // Test line intersection
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (x = minX; x <= maxX; x += grid->getDelta())
                {
                    if (mesh->intersectLine((x - grid->getDelta()), y, z, x, y, z)) 
                        break;
                    grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                }
            }
        }

        // Test line intersection from opposite direction
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (x = maxX; x >= minX; x -= grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        if (mesh->intersectLine((x + grid->getDelta()), y, z, x, y, z))
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (x = minX; x <= maxX; x += grid->getDelta())
            {
                for (y = minY; y <= maxY; y += grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        if (mesh->intersectLine(x, (y - grid->getDelta()), z, x, y, z)) 
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection from opposite direction
        for (z = minZ; z <= maxZ; z += grid->getDelta())
        {
            for (x = minX; x <= maxX; x += grid->getDelta())
            {
                for (y = maxY; y >= minY; y -= grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        if (mesh->intersectLine(x, (y + grid->getDelta()), z, x, y, z))
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection
        for (x = minX; x <= maxX; x += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (z = minZ; z <= maxZ; z += grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        if (mesh->intersectLine(x, y, (z - grid->getDelta()), x, y, z)) 
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        // Test line intersection from opposite direction
        for (x = minX; x <= maxX; x += grid->getDelta())
        {
            for (y = minY; y <= maxY; y += grid->getDelta())
            {
                for (z = maxZ; z >= minZ; z -= grid->getDelta())
                {
                    if (!grid->isNode(grid->transCoordToIndex(x, y, z), OuterType))
                    {
                        if (mesh->intersectLine(x, y, (z + grid->getDelta()), x, y, z)) 
                            break;
                        grid->setNodeTo(grid->transCoordToIndex(x, y, z), OuterType);
                    }
                }
            }
        }

        delete mesh;
}



void PointUnderTriangleStrategy::doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char innerType, char outerType)
{
#pragma omp parallel for
    for (long i = 0; i < triangularMesh->size; i++)
        this->meshReverse(triangularMesh->triangles[i], grid, innerType);

    this->findInsideNodes(grid, innerType);

#pragma omp parallel for
    for (int i = 0; i < (int)grid->getSize(); i++)
        this->setNegativeDirBorderTo(grid, i, innerType);
}

void PointUnderTriangleStrategy::meshReverse(Triangle& triangle, GridImp* grid, char innerType)
{
    auto box = grid->getBoundingBoxOnNodes(triangle);

    const real delta = grid->getDelta();
    triangle.initalLayerThickness(delta);

    for (real x = box.minX; x <= box.maxX; x += delta)
    {
        for (real y = box.minY; y <= box.maxY; y += delta)
        {
            for (real z = box.minZ; z <= box.maxZ; z += delta)
            {
                const uint index = grid->transCoordToIndex(x, y, z);

                const Vertex point(x, y, z);

                const char pointValue = triangle.isUnderFace(point);

                if (pointValue == vf::gpu::NEGATIVE_DIRECTION_BORDER)
                    grid->setNodeTo(index, vf::gpu::NEGATIVE_DIRECTION_BORDER);
                else if (pointValue == vf::gpu::INSIDE)
                    grid->setNodeTo(index, innerType);
            }
        }
    }
}

void PointUnderTriangleStrategy::findInsideNodes(GridImp* grid, char innerType)
{
    bool foundInsideNode = true;
    while (foundInsideNode)
    {
        foundInsideNode = false;
        for (uint index = 0; index < grid->getSize(); index++)
            this->setInsideNode(grid, index, foundInsideNode, innerType);
    }
}

void PointUnderTriangleStrategy::setInsideNode(GridImp* grid, const uint &index, bool &insideNodeFound, char innerType)
{
    if (grid->isNode(index, vf::gpu::NEGATIVE_DIRECTION_BORDER))
        return;

    if (!grid->isNode(index, innerType) && grid->nodeInNextCellIs(index, innerType))
    {
        grid->setNodeTo(index, innerType);
        insideNodeFound = true;
    }
}

void PointUnderTriangleStrategy::setNegativeDirBorderTo(GridImp* grid, const uint &index, char innerType)
{
    if (grid->isNode(index, vf::gpu::NEGATIVE_DIRECTION_BORDER))
        grid->setNodeTo(index, innerType);
}

//! \}
#if defined(VF_HAS_FAST_WINDING)
namespace
{
SPtr<GbTriFaceMesh3D> mergeWindingSurfaces(const std::vector<GbTriFaceMesh3D*>& surfaces)
{
    std::size_t totalNodes = 0;
    std::size_t totalTriangles = 0;
    std::vector<GbTriFaceMesh3D*> validSurfaces;
    validSurfaces.reserve(surfaces.size());

    for (const auto& surface : surfaces) {
        if (!surface)
            continue;

        auto* nodes = surface->getNodes();
        auto* tris  = surface->getTriangles();
        if (!nodes || !tris || nodes->empty() || tris->empty())
            continue;

        totalNodes += nodes->size();
        totalTriangles += tris->size();
        validSurfaces.push_back(surface);
    }

    if (validSurfaces.empty())
        return nullptr;

    auto* mergedNodes = new std::vector<GbTriFaceMesh3D::Vertex>();
    auto* mergedTriangles = new std::vector<GbTriFaceMesh3D::TriFace>();
    mergedNodes->reserve(totalNodes);
    mergedTriangles->reserve(totalTriangles);

    int nodeOffset = 0;
    for (auto* surface : validSurfaces) {
        auto* nodes = surface->getNodes();
        auto* tris  = surface->getTriangles();
        if (!nodes || !tris)
            continue;

        mergedNodes->insert(mergedNodes->end(), nodes->begin(), nodes->end());
        for (const auto& tri : *tris) {
            mergedTriangles->emplace_back(
                tri.getIndexVertex1() + nodeOffset,
                tri.getIndexVertex2() + nodeOffset,
                tri.getIndexVertex3() + nodeOffset);
        }

        nodeOffset += static_cast<int>(nodes->size());
    }

    if (mergedTriangles->empty()) {
        delete mergedNodes;
        delete mergedTriangles;
        return nullptr;
    }

    return std::make_shared<GbTriFaceMesh3D>(
        "fast_winding_merged_surface",
        mergedNodes,
        mergedTriangles,
        GbTriFaceMesh3D::KDTREE_SAHPLIT,
        false);
}
}

void FastWindingDiscretizationStrategy::doDiscretize(TriangularMesh* triangularMesh, GridImp* grid, char /*innerType*/, char /*outerType*/)
{
    if (!triangularMesh || !grid)
        return;

    triangularMesh->generateGbTriFaceMesh3D();
    auto surface = triangularMesh->VF_GbTriFaceMesh3D;
    if (!surface) {
        VF_LOG_WARNING("Fast winding strategy requested but no GbTriFaceMesh3D is available.");
        grid->setActiveWindingSurface(nullptr);
        return;
    }

    grid->setActiveWindingSurface(surface);

    vf::gpu::grid_winding::markSolidsWithFastWinding(*grid, *surface, accuracyScale, threshold, tolerance);
}

void FastWindingDiscretizationStrategy::appendFastWindingQSurfaces(
    GridImp* grid, std::vector<SPtr<GbTriFaceMesh3D>>& surfaces) const
{
    if (!grid)
        return;

    auto surface = grid->getActiveWindingSurface();
    if (surface) {
        VF_LOG_TRACE("FWN Qs: solid BC nodes = {}", grid->getNumberOfSolidBoundaryNodes());
        surfaces.push_back(surface);
    }
}

void FastWindingDiscretizationStrategy::computeFastWindingQs(
    GridImp* grid, const std::vector<SPtr<GbTriFaceMesh3D>>& surfaces) const
{
    if (!grid || surfaces.empty())
        return;

    // Keep only valid triangle surfaces for the final Q pass
    std::vector<GbTriFaceMesh3D*> surfacePtrs;
    surfacePtrs.reserve(surfaces.size());
    for (const auto& surface : surfaces) {
        if (!surface)
            continue;

        auto* nodes = surface->getNodes();
        auto* tris  = surface->getTriangles();
        if (!nodes || !tris || tris->empty())
            continue;

        surfacePtrs.push_back(surface.get());
    }

    if (surfacePtrs.empty())
        return;

    // Try to merge surfaces into one helper mesh to reduce kd tree work
    const std::size_t sourceSurfaceCount = surfacePtrs.size();
    auto mergedSurface = mergeWindingSurfaces(surfacePtrs);
    if (mergedSurface) {
        surfacePtrs.assign(1, mergedSurface.get());
        auto* mergedTris = mergedSurface->getTriangles();
        const std::size_t mergedTriangleCount = mergedTris ? mergedTris->size() : 0;
        VF_LOG_TRACE("FWN helper: merged {} surfaces into one kd-tree surface (triangles={})",
                     sourceSurfaceCount,
                     mergedTriangleCount);
    } else {
        VF_LOG_TRACE("FWN helper: using {} individual surfaces for Q computation", sourceSurfaceCount);
    }

    // Update BC_SOLID boundary nodes and rebuild qIndices before we compute Q values
    grid->rebuildBoundaryQIndices();
    // Allocate and reset Q storage for the current boundary node set
    grid->ensureQStorageAllocated();

    vf::gpu::grid_winding::BoundaryProcessingConfig boundaryConfig;
    vf::basics::Timer passTimer;
    passTimer.start();
    vf::grid_winding::SubgridDistanceStats stats;
    // Compute Q values from boundary links and triangle intersections using the helper
    vf::gpu::grid_winding::processSolidBoundaryLinksMultithreaded(
        *grid,
        surfacePtrs,
        &stats,
        boundaryConfig.parallelize,
        boundaryConfig.parallelThreshold,
        &surfacePtrs);
    const double passDuration = passTimer.getCurrentRuntimeInSeconds();
    VF_LOG_TRACE("FWN helper: finished Q computation (surfaces={}, hits={}, missing={}, time={} sec)",
                 surfacePtrs.size(),
                 stats.hitLinks,
                 stats.missingLinks.size(),
                 passDuration);

    // Fill remaining missing Q values with a fallback near solid neighbours
    grid->fillMissingQsAlongSolidNeighbours(static_cast<real>(vf::grid_winding::defaultMissingQ()));
    VF_LOG_TRACE("FWN Qs (post-pass): solid BC nodes = {}", grid->getNumberOfSolidBoundaryNodes());
}

}

#endif
