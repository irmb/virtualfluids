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
//! \author Hussein Alihussein
//! \note Generated with assistance from ChatGPT (GPT-5).
//=======================================================================================
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include <basics/geometry3d/winding/GridWindingFastDefaults.h>
#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/utilities/UbMath.h>

#if defined(VF_HAS_FAST_WINDING)
#include <igl/FastWindingNumberForSoups.h>
#endif

namespace vf::grid_winding
{
#if defined(VF_HAS_FAST_WINDING)
//! \brief Mark grid nodes as solid using the fast-winding approximation.
//!
//! \param accuracyScale Speed/accuracy control passed to libigl's `computeSolidAngle(...)`.
//!        Higher values are usually faster but less exact for small details.
//!        Lower values are usually slower but more exact.
//! \param threshold Normalized winding cutoff after scaling by `1 / (4*pi)`.
//!        `0.5` is a common inside/outside split for closed, well-oriented surfaces.
//!        Lower values mark more nodes as solid. Higher values mark fewer nodes as solid.
//! \param tolerance Comparison margin around the threshold to reduce noise near the surface.
//!        If `0` is passed, the shared default tolerance is used.
template <typename Accessor>
void markSolidsWithFastWinding(Accessor& accessor, GbTriFaceMesh3D& surface,
                               float accuracyScale = FastWindingDefaultAccuracyScale,
                               float threshold = FastWindingDefaultThreshold,
                               float tolerance = FastWindingDefaultTolerance)
{
    const auto* nodes = surface.getNodes();
    const auto* tris  = surface.getTriangles();
    if (!nodes || !tris || nodes->empty() || tris->empty())
        return;

    using FastVector    = igl::FastWindingNumber::HDK_Sample::UT_Vector3T<float>;
    using FastSolidAngle = igl::FastWindingNumber::HDK_Sample::UT_SolidAngle<float, float>;

    std::vector<FastVector> meshPositions;
    meshPositions.reserve(nodes->size());
    for (const auto& node : *nodes) {
        FastVector p;
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
                     meshPositions.data(), FastWindingDefaultExpansionOrder);

    // `0` means "use project default tolerance" to avoid disabling the safety margin by accident.
    const float     effectiveThreshold   = threshold;
    const float     effectiveTolerance   =
        ub_math::zero(tolerance) ? FastWindingDefaultTolerance : tolerance;
    const float     invFourPi            = 1.0f / (4.0f * static_cast<float>(std::acos(-1.0)));

    const std::array<double, 3> meshMin{ surface.getX1Minimum(), surface.getX2Minimum(), surface.getX3Minimum() };
    const std::array<double, 3> meshMax{ surface.getX1Maximum(), surface.getX2Maximum(), surface.getX3Maximum() };
    const double                faceEps = computeFaceEpsilon(meshMin, meshMax, accessor.cellSizeHint());

    accessor.forEachNode([&](typename Accessor::Node& node) {
        const auto& p = node.position;
        const double px = p.X1();
        const double py = p.X2();
        const double pz = p.X3();

        if (ub_math::less(px, meshMin[0] - faceEps) || ub_math::greater(px, meshMax[0] + faceEps) ||
            ub_math::less(py, meshMin[1] - faceEps) || ub_math::greater(py, meshMax[1] + faceEps) ||
            ub_math::less(pz, meshMin[2] - faceEps) || ub_math::greater(pz, meshMax[2] + faceEps))
            return;

        FastVector query;
        query[0] = static_cast<float>(px);
        query[1] = static_cast<float>(py);
        query[2] = static_cast<float>(pz);

        // `accuracyScale` is the libigl fast-winding accuracy/speed parameter from the cited method.
        const float winding = fastWinding.computeSolidAngle(query, accuracyScale) * invFourPi;
        if (ub_math::greaterEqual(winding, effectiveThreshold - effectiveTolerance))
            accessor.setSolid(node, winding);
    });
}

//! \brief Convenience overload for const surface meshes.
//!
//! Defaults are repeated here intentionally because this is a separate overload in C++.
template <class Accessor>
inline void markSolidsWithFastWinding(Accessor& accessor, const GbTriFaceMesh3D& surface,
                                      float accuracyScale = FastWindingDefaultAccuracyScale,
                                      float threshold = FastWindingDefaultThreshold,
                                      float tolerance = FastWindingDefaultTolerance)
{
    auto& nc = const_cast<GbTriFaceMesh3D&>(surface);
    markSolidsWithFastWinding(accessor, nc, accuracyScale, threshold, tolerance);
}
#endif

} // namespace vf::grid_winding
