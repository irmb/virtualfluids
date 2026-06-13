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
//! \brief Helper functions for winding-based solid marking on GPU triangular mesh grids.
//! \note Generated with the assistance of OpenAI Codex (GPT-5). Reviewed and adapted by author.
//=======================================================================================
#pragma once

#include <array>

#include <basics/geometry3d/winding/GridWindingFastDefaults.h>
#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>
#include <basics/geometry3d/winding/GridWindingClassification.h>
// #include <basics/geometry3d/winding/GeneralizedWindingNumber.h>

#include "grid/GridImp.h"
#include <grid/NodeValues.h>

namespace vf::gpu::grid_winding
{
// Adapter that exposes GridImp nodes to the shared fast-winding solid marker logic.
class GpuWindingSolidAccessor
{
public:
    using Vec3 = vf::grid_winding::Vec3;

    struct Node
    {
        Vec3 position;
        GridImp *grid{ nullptr };
        uint index{ 0 };
    };

    explicit GpuWindingSolidAccessor(GridImp &grid)
        : grid_(grid), cellSizeHint_(static_cast<double>(grid.getDelta()))
    {
    }

    template <typename Fn>
    void forEachNode(Fn &&fn)
    {
        const uint size = grid_.getSize();
        for (uint index = 0; index < size; ++index) {
            real x{ 0.0 };
            real y{ 0.0 };
            real z{ 0.0 };
            grid_.transIndexToCoords(index, x, y, z);

            Node node;
            node.position = Vec3{ static_cast<double>(x), static_cast<double>(y), static_cast<double>(z) };
            node.grid     = &grid_;
            node.index    = index;

            fn(node);
        }
    }

    void setSolid(const Node &node, float /*windingValue*/ = 0.0f) const
    {
        if (!node.grid)
            return;

        const auto current = node.grid->getFieldEntry(node.index);
        if (current == FLUID || current == BC_SOLID || current == STOPPER_SOLID)
            node.grid->setFieldEntry(node.index, INVALID_SOLID);
    }

    [[nodiscard]] double cellSizeHint() const
    {
        return cellSizeHint_;
    }

private:
    GridImp &grid_;
    double cellSizeHint_{ 0.0 };
};

#if defined(VF_HAS_FAST_WINDING)
//! \brief GPU wrapper around the shared fast-winding solid marker.
//!
//! \note Defaults are repeated here on purpose. In C++, default arguments are not
//! inherited from the called function, so each direct entry point needs its own defaults.
inline void markSolidsWithFastWinding(
    GridImp &grid, const GbTriFaceMesh3D &surface,
    float accuracyScale = vf::grid_winding::FastWindingDefaultAccuracyScale,
    float threshold = vf::grid_winding::FastWindingDefaultThreshold,
    float tolerance = vf::grid_winding::FastWindingDefaultTolerance)
{
    GpuWindingSolidAccessor accessor(grid);
    vf::grid_winding::markSolidsWithFastWinding(accessor, surface, accuracyScale, threshold, tolerance);
}
#endif

} // namespace vf::grid_winding::gpu
