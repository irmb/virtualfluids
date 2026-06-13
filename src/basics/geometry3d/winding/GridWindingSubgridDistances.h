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
//! \note Generated with assistance from ChatGPT (GPT-5).
//=======================================================================================
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <vector>

#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/geometry3d/GbVector3D.h>
#include <basics/utilities/UbKeys.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>

#include <geometry3d/KdTree/KdRay.h>
#include <geometry3d/KdTree/KdTree.h>

namespace vf::grid_winding
{
using Vec3 = GbVector3D;
using TriFace = GbTriFaceMesh3D::TriFace;
using TriVertex = GbTriFaceMesh3D::Vertex;
using TriangleSurfaces = std::vector<GbTriFaceMesh3D *>;

inline double defaultMissingQ()
{
    return 0.5;
}

inline double computeFaceEpsilon(const std::array<double, 3>& meshMin, const std::array<double, 3>& meshMax,
                                 double                       cellScale = 0.0)
{
    double maxAbs = 1.0;
    for (double value : meshMin)
        maxAbs = std::max(maxAbs, std::abs(value));
    for (double value : meshMax)
        maxAbs = std::max(maxAbs, std::abs(value));

    constexpr double scale   = 128.0;
    double           epsilon = maxAbs * std::numeric_limits<double>::epsilon() * scale;
    if (ub_math::greater(cellScale, 0.0))
        epsilon = std::max(epsilon, cellScale * 1.0e-6);
    return epsilon;
}

inline void accumulateBounds(Vec3& minV, Vec3& maxV, const Vec3& p)
{
    minV.X1() = std::min(minV.X1(), p.X1());
    minV.X2() = std::min(minV.X2(), p.X2());
    minV.X3() = std::min(minV.X3(), p.X3());
    maxV.X1() = std::max(maxV.X1(), p.X1());
    maxV.X2() = std::max(maxV.X2(), p.X2());
    maxV.X3() = std::max(maxV.X3(), p.X3());
}

inline void segmentBounds(const Vec3& a, const Vec3& b, Vec3& outMin, Vec3& outMax)
{
    outMin.X1() = std::min(a.X1(), b.X1());
    outMin.X2() = std::min(a.X2(), b.X2());
    outMin.X3() = std::min(a.X3(), b.X3());
    outMax.X1() = std::max(a.X1(), b.X1());
    outMax.X2() = std::max(a.X2(), b.X2());
    outMax.X3() = std::max(a.X3(), b.X3());
}

inline bool overlapsSegmentAabb(const Vec3& aMin, const Vec3& aMax, const Vec3& bMin, const Vec3& bMax,
                                double      eps = 0.0)
{
    return (ub_math::lessEqual(aMin.X1(), bMax.X1() + eps) && ub_math::greaterEqual(aMax.X1() + eps, bMin.X1())) &&
           (ub_math::lessEqual(aMin.X2(), bMax.X2() + eps) && ub_math::greaterEqual(aMax.X2() + eps, bMin.X2())) &&
           (ub_math::lessEqual(aMin.X3(), bMax.X3() + eps) && ub_math::greaterEqual(aMax.X3() + eps, bMin.X3()));
}


struct SubgridDistanceStats
{
    struct MissingLink
    {
        std::array<double, 3> origin{};
        std::array<double, 3> neighbour{};
        int                   direction{ 0 };
        float                 windingNumber{ std::numeric_limits<float>::quiet_NaN() };
        float                 endpointKind{ 0.0f };
    };

    std::size_t              totalLinks{ 0 };
    std::size_t              hitLinks{ 0 };
    std::vector<MissingLink> missingLinks;
    std::size_t              promotedFluidNodes{ 0 };
    std::size_t              promotedZeroQNodes{ 0 };
    std::size_t              skippedNonLocalLinks{ 0 };
};

//==============================================================================
// Diagnostics data containers (common CPU/GPU)
//==============================================================================

struct MissingLinkCollection
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2>   lines;
    std::vector<float>         directions;
    std::vector<float>         windingNumbers;
    std::vector<float>         endpointKind;
    std::size_t                total{ 0 };

    [[nodiscard]] bool empty() const { return total == 0; }
};

struct QLineCollection
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2>   lines;
    std::vector<float>         directions;
    std::vector<float>         qValues;
    std::vector<UbTupleFloat3> surfacePoints;
    std::size_t                total{ 0 };

    [[nodiscard]] bool empty() const { return total == 0; }
};

namespace cpu
{
struct BoundaryProcessingConfig
{
    bool   writeMissingLinks{ true };
    bool   writeQLines{ true };
    // Fast-winding solid classification threshold (CPU).
    float  fastWindingThreshold{ 0.5f };
};

struct BoundaryProcessingResult
{
    vf::grid_winding::SubgridDistanceStats  stats;
    vf::grid_winding::MissingLinkCollection missingLinks;
    vf::grid_winding::QLineCollection       qLines;
};
} // namespace cpu

namespace detail
{
constexpr double kEps = 1.0e-12;
constexpr double kBaryTol = 1.0e-6;

struct TriangleInfo
{
    Vec3   v0;
    Vec3   e1;
    Vec3   e2;
    Vec3   normal;
    double planeD{ 0.0 };
    Vec3   bboxMin;
    Vec3   bboxMax;
};

inline TriangleInfo buildTriangleInfo(const TriFace& tri, const std::vector<TriVertex>& nodes)
{
    TriangleInfo info;

    const auto& n0 = nodes[tri.getIndexVertex1()];
    const auto& n1 = nodes[tri.getIndexVertex2()];
    const auto& n2 = nodes[tri.getIndexVertex3()];

    info.v0 = Vec3(static_cast<double>(n0.x), static_cast<double>(n0.y), static_cast<double>(n0.z));
    Vec3 v1(static_cast<double>(n1.x), static_cast<double>(n1.y), static_cast<double>(n1.z));
    Vec3 v2(static_cast<double>(n2.x), static_cast<double>(n2.y), static_cast<double>(n2.z));

    info.e1     = v1 - info.v0;
    info.e2     = v2 - info.v0;
    info.normal = info.e1.Cross(info.e2);
    info.normal.Normalize();
    info.planeD = info.normal.Dot(info.v0);

    info.bboxMin = info.bboxMax = info.v0;
    accumulateBounds(info.bboxMin, info.bboxMax, v1);
    accumulateBounds(info.bboxMin, info.bboxMax, v2);

    return info;
}

inline bool pointInTriangle(const TriangleInfo& triangle, const Vec3& p, double baryTol = kBaryTol)
{
    const Vec3 v0 = triangle.e1;
    const Vec3 v1 = triangle.e2;
    const Vec3 v2 = p - triangle.v0;

    const double dot00 = v0.Dot(v0);
    const double dot01 = v0.Dot(v1);
    const double dot02 = v0.Dot(v2);
    const double dot11 = v1.Dot(v1);
    const double dot12 = v1.Dot(v2);

    const double denom = dot00 * dot11 - dot01 * dot01;
    if (ub_math::less(std::fabs(denom), kEps))
        return false;

    const double invDenom = 1.0 / denom;
    const double u        = (dot11 * dot02 - dot01 * dot12) * invDenom;
    const double v        = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return ub_math::greaterEqual(u, -baryTol) && ub_math::greaterEqual(v, -baryTol) &&
           ub_math::lessEqual(u + v, 1.0 + baryTol);
}

inline bool pointOnTriangle(const TriangleInfo& triangle, const Vec3& p, double planeTol, double baryTol = 1.0e-8)
{
    const double planeDist = triangle.normal.Dot(p) - triangle.planeD;
    if (ub_math::greater(std::fabs(planeDist), planeTol))
        return false;
    return pointInTriangle(triangle, p, baryTol);
}

inline bool pointNearBounds(const TriangleInfo& triangle, const Vec3& p, double tol)
{
    return ub_math::greaterEqual(p.X1(), triangle.bboxMin.X1() - tol) &&
           ub_math::lessEqual(p.X1(), triangle.bboxMax.X1() + tol) &&
           ub_math::greaterEqual(p.X2(), triangle.bboxMin.X2() - tol) &&
           ub_math::lessEqual(p.X2(), triangle.bboxMax.X2() + tol) &&
           ub_math::greaterEqual(p.X3(), triangle.bboxMin.X3() - tol) &&
           ub_math::lessEqual(p.X3(), triangle.bboxMax.X3() + tol);
}

inline int endpointSurfaceKind(const TriangleSurfaces& surfaces, const Vec3& origin, const Vec3& neighbour,
                               double planeTol, double baryTol = 1.0e-8)
{
    if (surfaces.empty())
        return 0;

    int kind = 0;
    const double bboxTol = planeTol;

    for (auto* surface : surfaces) {
        if (!surface)
            continue;
        auto* nodes = surface->getNodes();
        auto* tris = surface->getTriangles();
        if (!nodes || !tris)
            continue;

        for (const auto& tri : *tris) {
            TriangleInfo info = buildTriangleInfo(tri, *nodes);

            if (!(kind & 0x1) && pointNearBounds(info, origin, bboxTol) &&
                pointOnTriangle(info, origin, planeTol, baryTol)) {
                kind |= 0x1;
            }
            if (!(kind & 0x2) && pointNearBounds(info, neighbour, bboxTol) &&
                pointOnTriangle(info, neighbour, planeTol, baryTol)) {
                kind |= 0x2;
            }
            if (kind == 3)
                return kind;
        }
    }

    return kind;
}

inline int endpointSurfaceKind(const std::vector<TriangleInfo>& triangles, const Vec3& origin, const Vec3& neighbour,
                               double planeTol, double baryTol = 1.0e-8)
{
    if (triangles.empty())
        return 0;

    int kind = 0;
    const double bboxTol = planeTol;

    for (const auto& tri : triangles) {
        if (!(kind & 0x1) && pointNearBounds(tri, origin, bboxTol) &&
            pointOnTriangle(tri, origin, planeTol, baryTol)) {
            kind |= 0x1;
        }
        if (!(kind & 0x2) && pointNearBounds(tri, neighbour, bboxTol) &&
            pointOnTriangle(tri, neighbour, planeTol, baryTol)) {
            kind |= 0x2;
        }
        if (kind == 3)
            return kind;
    }

    return kind;
}

inline bool intersectTriangle(const TriangleInfo& tri, const Vec3& origin, const Vec3& dir, double maxT,
                              double& outT)
{
    const Vec3  pvec{ dir.X2() * tri.e2.X3() - dir.X3() * tri.e2.X2(),
                      dir.X3() * tri.e2.X1() - dir.X1() * tri.e2.X3(),
                      dir.X1() * tri.e2.X2() - dir.X2() * tri.e2.X1() };
    const double det        = tri.e1.Dot(pvec);
    const double segmentTol = std::max(1.0e-10, 1.0e-6 * std::max(1.0, maxT));
    const double uvTol      = kBaryTol;

    if (ub_math::less(std::fabs(det), kEps)) {
        const Vec3  endPoint{ origin.X1() + dir.X1() * maxT,
                              origin.X2() + dir.X2() * maxT,
                              origin.X3() + dir.X3() * maxT };
        const double originPlane = tri.normal.Dot(origin) - tri.planeD;
        const double endPlane    = tri.normal.Dot(endPoint) - tri.planeD;

        const double planeTol = segmentTol;

        if (ub_math::lessEqual(std::fabs(originPlane), planeTol) && pointInTriangle(tri, origin)) {
            outT = 0.0;
            return true;
        }

        if (ub_math::lessEqual(std::fabs(endPlane), planeTol) && pointInTriangle(tri, endPoint)) {
            outT = maxT;
            return true;
        }

        return false;
    }

    const double invDet = 1.0 / det;

    const Vec3  tvec{ origin.X1() - tri.v0.X1(),
                      origin.X2() - tri.v0.X2(),
                      origin.X3() - tri.v0.X3() };
    const double u = tvec.Dot(pvec) * invDet;
    if (ub_math::less(u, -uvTol) || ub_math::greater(u, 1.0 + uvTol))
        return false;

    const Vec3  qvec{ tvec.X2() * tri.e1.X3() - tvec.X3() * tri.e1.X2(),
                      tvec.X3() * tri.e1.X1() - tvec.X1() * tri.e1.X3(),
                      tvec.X1() * tri.e1.X2() - tvec.X2() * tri.e1.X1() };
    const double v = dir.Dot(qvec) * invDet;
    if (ub_math::less(v, -uvTol) || ub_math::greater(u + v, 1.0 + uvTol))
        return false;

    const double t = tri.e2.Dot(qvec) * invDet;
    if (ub_math::less(t, -segmentTol) || ub_math::greater(t, maxT + segmentTol))
        return false;

    if (ub_math::lessEqual(t, segmentTol)) {
        outT = 0.0;
        return true;
    }

    if (ub_math::greaterEqual(t, maxT - segmentTol)) {
        outT = maxT;
        return true;
    }

    outT = t;
    return true;
}

} // namespace detail

inline bool pointOnSurface(GbTriFaceMesh3D &surface, const Vec3 &point, double eps)
{
    return surface.intersectLine(point.X1() - eps, point.X2(), point.X3(),
                                 point.X1() + eps, point.X2(), point.X3()) ||
           surface.intersectLine(point.X1(), point.X2() - eps, point.X3(),
                                 point.X1(), point.X2() + eps, point.X3()) ||
           surface.intersectLine(point.X1(), point.X2(), point.X3() - eps,
                                 point.X1(), point.X2(), point.X3() + eps);
}

inline int endpointSurfaceKind(GbTriFaceMesh3D &surface, const Vec3 &origin, const Vec3 &neighbour, double eps)
{
    int kind = 0;
    if (pointOnSurface(surface, origin, eps))
        kind |= 0x1;
    if (pointOnSurface(surface, neighbour, eps))
        kind |= 0x2;
    return kind;
}

template <typename PromoteSolidFn, typename ExistingQFn, typename SetQFn, typename MarkSurfaceFn,
          typename GetWindingFn, typename HandleMissingFn>
inline bool processBoundaryLink(const TriangleSurfaces& surfaces, const Vec3& origin, double dx, int cx,
                                int cy, int cz, int direction, PromoteSolidFn&& promoteSolid,
                                ExistingQFn&& existingQ, SetQFn&& setQ, MarkSurfaceFn&& markSurface,
                                GetWindingFn&& getWindingNumber, HandleMissingFn&& handleMissing,
                                SubgridDistanceStats* stats)
{
    Vec3 neighbour{ origin.X1() + dx * static_cast<double>(cx),
                    origin.X2() + dx * static_cast<double>(cy),
                    origin.X3() + dx * static_cast<double>(cz) };

    Vec3 segMin;
    Vec3 segMax;
    segmentBounds(origin, neighbour, segMin, segMax);

    Vec3        dirVec      = neighbour - origin;
    const double maxDistance = dirVec.Length();
    if (ub_math::less(maxDistance, detail::kEps))
        return false;

    dirVec.Normalize();

    double bestT = std::numeric_limits<double>::max();
    int    bestPatchIndex = -1;
    bool   hasHit         = false;

    for (auto* surface : surfaces) {
        if (!surface)
            continue;
        auto* nodes = surface->getNodes();
        auto* tris = surface->getTriangles();
        if (!nodes || !tris)
            continue;

        for (const auto& tri : *tris) {
            detail::TriangleInfo info = detail::buildTriangleInfo(tri, *nodes);

            if (!overlapsSegmentAabb(segMin, segMax, info.bboxMin, info.bboxMax))
                continue;

            const double signedDistance = info.normal.Dot(origin) - info.planeD;
            const double bandLimit = 2.1 * maxDistance;
            if (ub_math::less(signedDistance, -bandLimit))
                continue;

            double hitT = 0.0;
            if (detail::intersectTriangle(info, origin, dirVec, maxDistance, hitT)) {
                hasHit = true;
                if (ub_math::less(hitT, bestT))
                    bestT = hitT;
            }
        }
    }

    if (!hasHit) {
        typename SubgridDistanceStats::MissingLink missing;
        missing.origin    = { origin.X1(), origin.X2(), origin.X3() };
        missing.neighbour = { neighbour.X1(), neighbour.X2(), neighbour.X3() };
        missing.direction = direction;
        const bool handled = handleMissing(direction, origin, neighbour, missing);
        if (!handled && promoteSolid() && stats)
            ++stats->promotedFluidNodes;
        if (stats) {
            auto windingInfo = getWindingNumber();
            if (windingInfo.first)
                missing.windingNumber = static_cast<float>(windingInfo.second);
            stats->missingLinks.push_back(missing);
        }
        return false;
    }

    if (stats)
        ++stats->hitLinks;

    if (ub_math::equal(bestT, std::numeric_limits<double>::max()))
        return true;

    const double rawQ          = bestT / maxDistance;
    constexpr double qZeroTol    = 5.0e-6;
    constexpr double qUpperTol   = 5.0e-6;
    constexpr double qSurfaceTol = 1.0e-6;

    if (ub_math::less(rawQ, -qZeroTol))
        return true;

    double q = rawQ;

    if (ub_math::lessEqual(std::fabs(rawQ), qZeroTol)) {
        q = 0.0;
    } else if (ub_math::greater(rawQ, 1.0) && ub_math::lessEqual(rawQ, 1.0 + qUpperTol)) {
        q = 1.0;
    }

    q = std::clamp(q, 0.0, 1.0);

    if (ub_math::lessEqual(std::fabs(q - 1.0), qSurfaceTol))
        markSurface(q, bestPatchIndex);

    const double existing = existingQ();
    if (ub_math::greater(existing, 0.0) && ub_math::lessEqual(existing, q))
        return true;

    setQ(q, bestPatchIndex);
    return true;
}

struct BoundaryLink
{
    Vec3 origin;
    Vec3 neighbour;
    int  direction{ 0 };
};

struct QComputationConfig
{
    double qZeroTol{ 5.0e-6 };
    double qUpperTol{ 5.0e-6 };
    double qSurfaceTol{ 1.0e-6 };
};

//! \brief Ray-intersection handler that returns the closest hit factor along a ray.
class KdClosestIntersectionHandler : public kd_tree::RayIntersectionHandler<double>
{
public:
    explicit KdClosestIntersectionHandler(double maxFactor)
        : maxFactor_(maxFactor)
        , bestFactor_(maxFactor)
        , hit_(false)
    {
    }

    void reset(double maxFactor)
    {
        maxFactor_  = maxFactor;
        bestFactor_ = maxFactor;
        hit_        = false;
    }

    [[nodiscard]] bool hit() const { return hit_; }

    [[nodiscard]] double bestFactor() const { return bestFactor_; }

    int intersectRay(const kd_tree::Ray<double> &ray, kd_tree::Node<double> &parent, kd_tree::Node<double> *&child1,
                     kd_tree::Node<double> *&child2, std::set<ub_keys::Key3<int>> &mailbox) const override
    {
        const int bboxState = parent.intersectRayBoundingBox(ray);
        if (bboxState != kd_tree::Intersection::INTERSECTION)
            return kd_tree::Intersection::NO_INTERSECTION;

        if (parent.isLeaf()) {
            const auto &triFacesPtr = parent.getTriFaces();
            if (!triFacesPtr)
                return kd_tree::Intersection::NO_INTERSECTION;

            auto       &triFaces = *triFacesPtr;
            auto       &nodes    = parent.getNodes();

            const Vec3 origin(ray.originX, ray.originY, ray.originZ);
            const Vec3 dir(ray.directionX, ray.directionY, ray.directionZ);

            for (auto &triFace : triFaces) {
                ub_keys::Key3<int> key(triFace.getIndexVertex1(), triFace.getIndexVertex2(),
                                       triFace.getIndexVertex3());
                if (mailbox.find(key) != mailbox.end())
                    continue;
                mailbox.insert(key);

                double hitT = 0.0;
                const auto info = detail::buildTriangleInfo(triFace, nodes);
                if (!detail::intersectTriangle(info, origin, dir, maxFactor_, hitT))
                    continue;

                if (ub_math::zero(hitT)) {
                    hit_        = true;
                    bestFactor_ = 0.0;
                    return kd_tree::Intersection::ON_BOUNDARY;
                }

                if (ub_math::less(hitT, bestFactor_)) {
                    bestFactor_ = hitT;
                    hit_        = true;
                }
            }

            return hit_ ? kd_tree::Intersection::INTERSECTION : kd_tree::Intersection::NO_INTERSECTION;
        }

        int result = kd_tree::Intersection::NO_INTERSECTION;
        if (child1) {
            const int childRes = child1->intersectRay(ray, *this, mailbox);
            if (childRes == kd_tree::Intersection::ON_BOUNDARY)
                return childRes;
            if (childRes == kd_tree::Intersection::INTERSECTION)
                result = kd_tree::Intersection::INTERSECTION;
        }
        if (child2) {
            const int childRes = child2->intersectRay(ray, *this, mailbox);
            if (childRes == kd_tree::Intersection::ON_BOUNDARY)
                return childRes;
            if (childRes == kd_tree::Intersection::INTERSECTION)
                result = kd_tree::Intersection::INTERSECTION;
        }

        return result;
    }

private:
    double              maxFactor_;
    mutable double      bestFactor_;
    mutable bool        hit_;
};

inline std::vector<kd_tree::Tree<double> *>
ensureKdTrees(const std::vector<GbTriFaceMesh3D *> &meshes)
{
    std::vector<kd_tree::Tree<double> *> kdTrees;
    kdTrees.reserve(meshes.size());
    for (auto *mesh : meshes) {
        if (!mesh)
            continue;
        if (auto *tree = mesh->ensureKdTree())
            kdTrees.push_back(tree);
    }
    return kdTrees;
}

inline std::vector<kd_tree::Tree<double> *>
ensureKdTrees(const std::vector<SPtr<GbTriFaceMesh3D>> &surfaces)
{
    std::vector<GbTriFaceMesh3D *> raw;
    raw.reserve(surfaces.size());
    for (const auto &surface : surfaces) {
        if (surface)
            raw.push_back(surface.get());
    }
    return ensureKdTrees(raw);
}

//! \brief Compute Q for a single boundary link using kd-tree based ray casting.
//!
//! Adapter must provide the following interface:
//!   bool promoteSolid(const BoundaryLink &link) const;
//!   double existingQ(const BoundaryLink &link) const;
//!   void setQ(const BoundaryLink &link, double q, int patchIndex) const;
//!   void markSurface(const BoundaryLink &link, double q, int patchIndex) const;
//!   std::pair<bool,double> windingNumber(const BoundaryLink &link) const;
//!   bool handleMissing(const BoundaryLink &link, SubgridDistanceStats::MissingLink &out) const;
template <class Adapter>
inline bool processBoundaryLinkKdTree(const std::vector<kd_tree::Tree<double> *> &kdTrees, Adapter &adapter,
                                      const BoundaryLink &link, SubgridDistanceStats *stats = nullptr,
                                      const QComputationConfig &cfg = {})
{
    using vf::grid_winding::detail::kEps;

    if (kdTrees.empty())
        return false;

    const Vec3 origin    = link.origin;
    const Vec3 neighbour = link.neighbour;

    Vec3        dirVec      = neighbour - origin;
    const double maxDistance = dirVec.Length();
    if (ub_math::less(maxDistance, kEps))
        return false;

    dirVec /= maxDistance;

    double bestT = std::numeric_limits<double>::max();
    bool   hit   = false;

    for (auto *tree : kdTrees) {
        if (!tree)
            continue;

        kd_tree::Ray<double> ray(origin.X1(), origin.X2(), origin.X3(),
                                 dirVec.X1(), dirVec.X2(), dirVec.X3());
        KdClosestIntersectionHandler handler(maxDistance);

        const int result = tree->intersectRay(ray, handler);
        if (result == kd_tree::Intersection::ON_BOUNDARY) {
            bestT = 0.0;
            hit   = true;
            break;
        }

        if (handler.hit() && ub_math::less(handler.bestFactor(), bestT)) {
            bestT = handler.bestFactor();
            hit   = true;
        }
    }

    if (!hit) {
        if (stats) {
            SubgridDistanceStats::MissingLink missing{};
            missing.origin    = { origin.X1(), origin.X2(), origin.X3() };
            missing.neighbour = { neighbour.X1(), neighbour.X2(), neighbour.X3() };
            missing.direction = link.direction;

            const bool handled = adapter.handleMissing(link, missing);
            if (!handled && adapter.promoteSolid(link))
                ++stats->promotedFluidNodes;

            const auto windingInfo = adapter.windingNumber(link);
            if (windingInfo.first)
                missing.windingNumber = static_cast<float>(windingInfo.second);

            stats->missingLinks.push_back(missing);
        } else {
            if (adapter.promoteSolid(link) && stats)
                ++stats->promotedFluidNodes;
        }

        return false;
    }

    if (stats)
        ++stats->hitLinks;

    if (ub_math::equal(bestT, std::numeric_limits<double>::max()))
        return true;

    const double rawQ = bestT / maxDistance;

    if (ub_math::less(rawQ, -cfg.qZeroTol))
        return true;

    double q = rawQ;

    if (ub_math::lessEqual(std::fabs(rawQ), cfg.qZeroTol)) {
        q = 0.0;
    } else if (ub_math::greater(rawQ, 1.0) && ub_math::lessEqual(rawQ, 1.0 + cfg.qUpperTol)) {
        q = 1.0;
    }

    q = std::clamp(q, 0.0, 1.0);

    if (ub_math::lessEqual(std::fabs(q - 1.0), cfg.qSurfaceTol))
        adapter.markSurface(link, q, -1);

    const double existing = adapter.existingQ(link);
    if (ub_math::greater(existing, 0.0) && ub_math::lessEqual(existing, q))
        return true;

    adapter.setQ(link, q, -1);
    return true;
}
} // namespace vf::grid_winding
