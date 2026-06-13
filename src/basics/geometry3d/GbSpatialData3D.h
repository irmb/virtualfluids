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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Henry Korb
//=======================================================================================
#ifndef GBPOINTCLOUD_H
#define GBPOINTCLOUD_H

#include "DataTypes.h"
#include "UbMath.h"
#include "logger/Logger.h"
#include <basics/utilities/UbObserver.h>

#include <GbObject3D.h>
#include <UbTuple.h>

#include <PointerDefinitions.h>
#include <cmath>
#include <constants/NumericConstants.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>

//////////////////////////////////////////////////////////////////////////
//!
//! \class GbPointCloud3D
//! \brief This class provides an object to interpolate spatially distributed data to new points.
//!
//////////////////////////////////////////////////////////////////////////
template <typename T>
class GbSpatialData3D
{
public:
    virtual T interpolateToPoint(real3 point) const = 0;
    virtual ~GbSpatialData3D() {};
};

//////////////////////////////////////////////////////////////////////////
//!
//! \class GbStructuredMesh3D
//! \brief This class provides an object to interpolate data on a regular structured mesh to new points.
template <typename T>
class GbStructuredMesh3D : public GbSpatialData3D<T>
{
public:
    enum class InterpolationStrategy { Trilinear };
    GbStructuredMesh3D(const real3& spacing, const real3& origin, const std::array<uint, 3>& nPoints,
                       const std::vector<T>& data, InterpolationStrategy interpolationStrategy)
        : spacing(spacing), origin(origin), nPoints(nPoints), data(data), interpolationStrategy(interpolationStrategy)
    {
        if (nPoints[0] * nPoints[1] * nPoints[2] != data.size())
            throw std::runtime_error("number of values does not match number of points.");
    }
    T interpolateToPoint(real3 point) const override
    {
        return applyTrilinearInterpolation(point);
    };
    ~GbStructuredMesh3D() override {};

private:
    T applyTrilinearInterpolation(const real3& point) const
    {
        if (isOutOfBounds(point))
            throw std::runtime_error("Point is out of bounds!");

        const auto [ix, iy, iz] = getIndexMMM(point);
        const real3 distances = (point - getCoordinates(ix, iy, iz)) / spacing;

        const bool endOfX = ix == nPoints[0] - 1;
        const bool endOfY = iy == nPoints[1] - 1;
        const bool endOfZ = iz == nPoints[2] - 1;

        if(endOfX && endOfY && endOfZ)
            return data.back();
        if (endOfX && endOfY)
            return ub_math::lerp(data[getIndex(ix, iy, iz)], data[getIndex(ix, iy, iz + 1)], distances.z);

        if (endOfX && endOfZ)
            return ub_math::lerp(data[getIndex(ix, iy, iz)], data[getIndex(ix, iy + 1, iz)], distances.y);

        if (endOfY && endOfZ)
            return ub_math::lerp(data[getIndex(ix, iy, iz)], data[getIndex(ix + 1, iy, iz)], distances.x);

        if (endOfX)
            return ub_math::lerp2(data[getIndex(ix, iy, iz)], data[getIndex(ix, iy + 1, iz)], data[getIndex(ix, iy, iz + 1)],
                                  data[getIndex(ix, iy + 1, iz + 1)], distances.y, distances.z);
        if (endOfY)
            return ub_math::lerp2(data[getIndex(ix, iy, iz)], data[getIndex(ix + 1, iy, iz)], data[getIndex(ix, iy, iz + 1)],
                                  data[getIndex(ix + 1, iy, iz + 1)], distances.x, distances.z);
        if (endOfZ)
            return ub_math::lerp2(data[getIndex(ix, iy, iz)], data[getIndex(ix + 1, iy + 1, iz)],
                                  data[getIndex(ix, iy + 1, iz)], data[getIndex(ix + 1, iy + 1, iz)], distances.x,
                                  distances.y);

        return ub_math::lerp3(data[getIndex(ix, iy, iz)], data[getIndex(ix + 1, iy, iz)], data[getIndex(ix, iy + 1, iz)],
                              data[getIndex(ix + 1, iy + 1, iz)], data[getIndex(ix, iy, iz + 1)],
                              data[getIndex(ix + 1, iy, iz + 1)], data[getIndex(ix, iy + 1, iz + 1)],
                              data[getIndex(ix + 1, iy + 1, iz + 1)], distances.x, distances.y, distances.z);
    }
    bool isBelowLowerX(real x) const
    {
        if (spacing.x == vf::basics::constant::c0o1 || nPoints[0] == 1)
            return false;
        return x < origin.x;
    }
    bool isBelowLowerY(real y) const
    {
        if (spacing.y == vf::basics::constant::c0o1 || nPoints[1] == 1)
            return false;
        return y < origin.y;
    }
    bool isBelowLowerZ(real z) const
    {
        if (spacing.z == vf::basics::constant::c0o1 || nPoints[2] == 1)
            return false;
        return z < origin.z;
    }

    bool isAboveUpperX(real x) const
    {
        if (spacing.x == vf::basics::constant::c0o1 || nPoints[0] == 1)
            return false;
        return x > origin.x + spacing.x * nPoints[0];
    }
    bool isAboveUpperY(real y) const
    {
        if (spacing.y == vf::basics::constant::c0o1 || nPoints[1] == 1)
            return false;
        return y > origin.y + spacing.y * nPoints[1];
    }
    bool isAboveUpperZ(real z) const
    {
        if (spacing.z == vf::basics::constant::c0o1 || nPoints[2] == 1)
            return false;
        return z > origin.z + spacing.z * nPoints[2];
    }

    bool isOutOfBounds(const real3& point) const
    {
        return isBelowLowerX(point.x) || isBelowLowerY(point.y) || isBelowLowerZ(point.z) || isAboveUpperX(point.x) ||
               isAboveUpperY(point.y) || isAboveUpperZ(point.z);
    }
    uint getNextLowerPoint(const real coord, const real origin, const real spacing, const uint numberOfPoints) const
    {
        if (spacing == vf::basics::constant::c0o1 || numberOfPoints == 1)
            return 0;
        if (coord < origin)
            return UINT32_MAX;
        return uint((coord - origin) / spacing);
    }
    uint getNextHigherPoint(const real coord, const real origin, const real spacing, const uint numberOfPoints) const
    {
        if (spacing == vf::basics::constant::c0o1 || numberOfPoints == 1)
            return 0;
        const uint index = (coord - origin) / spacing + 1;
        if (index > numberOfPoints)
            return UINT32_MAX;
        return index;
    }
    uint getNextLowerPointX(const real3& point) const
    {
        return getNextLowerPoint(point.x, origin.x, spacing.x, nPoints[0]);
    }
    uint getNextLowerPointY(const real3& point) const
    {
        return getNextLowerPoint(point.y, origin.y, spacing.y, nPoints[1]);
    }
    uint getNextLowerPointZ(const real3& point) const
    {
        return getNextLowerPoint(point.z, origin.z, spacing.z, nPoints[2]);
    }
    uint getNextHigherPointX(const real3& point) const
    {
        return getNextHigherPoint(point.x, origin.x, spacing.x, nPoints[0]);
    }
    uint getNextHigherPointY(const real3& point) const
    {
        return getNextHigherPoint(point.y, origin.y, spacing.y, nPoints[1]);
    }
    uint getNextHigherPointZ(const real3& point) const
    {
        return getNextHigherPoint(point.z, origin.z, spacing.z, nPoints[2]);
    }
    uint getIndex(const uint indexX, const uint indexY, const uint indexZ) const
    {
        return indexX + nPoints[0] * (indexY + nPoints[1] * indexZ);
    }
    std::array<uint, 3> getIndexMMM(const real3& point) const
    {
        return { getNextLowerPointX(point), getNextLowerPointY(point), getNextLowerPointZ(point) };
    }

    real getCoordinateX(const uint index) const
    {
        return origin.x + index * spacing.x;
    }
    real getCoordinateY(const uint index) const
    {
        return origin.y + index * spacing.y;
    }
    real getCoordinateZ(const uint index) const
    {
        return origin.z + index * spacing.z;
    }
    real3 getCoordinates(const uint indexX, const uint indexY, const uint indexZ) const
    {
        return { getCoordinateX(indexX), getCoordinateY(indexY), getCoordinateZ(indexZ) };
    }
    const real3 spacing, origin;
    const std::array<uint, 3> nPoints;
    std::vector<T> data;
    const InterpolationStrategy interpolationStrategy;
};

struct Node
{
    uint index;
    std::unique_ptr<Node> left, right;
    Node(uint index, std::unique_ptr<Node> left, std::unique_ptr<Node> right)
        : index(index), left(std::move(left)), right(std::move(right))
    {
    }
};

std::unique_ptr<Node> insertNodes(const std::vector<real3>& coordinates);
std::vector<uint> findNearestPoints(const Node* root, real3 point, uint nPoints, const std::vector<real3>& coordinates);
void applyToTree(Node* root, std::function<void(uint index, uint depth)> func);

//! \class GbPointCloud3D
//! \brief This class provides an object to interpolate irregularly distributed data to new points.
template <typename T>
class GbPointCloud3D : public GbSpatialData3D<T>
{
public:
    class InterpolationStrategy;
    class InverseDistanceWeighing;
    class NearestNeighbor;
    GbPointCloud3D(std::vector<real3> points, std::vector<T> data,
                   std::unique_ptr<InterpolationStrategy> interpolationStrategy, bool printTree = false)
        : coordinates(std::move(points)), values(std::move(data)), interpolationStrategy(std::move(interpolationStrategy))
    {
        if (points.size() != data.size())
            throw std::runtime_error("Size of points and data does not match!");
        fillTree();
        this->interpolationStrategy->setData(this);
        if (printTree)
            this->printTree();
    };
    ~GbPointCloud3D() override {};

    T interpolateToPoint(real3 point) const override
    {
        return this->interpolationStrategy->interpolateToPoint(point);
    };
    real3 getCoordinate(uint index) const
    {
        return this->coordinates[index];
    }
    T getValue(uint index) const
    {
        return this->values[index];
    }

    std::vector<uint> findNearestPoints(real3 point, uint nPoints) const
    {
        return ::findNearestPoints(root.get(), point, nPoints, coordinates);
    };
    uint findNearestPoint(real3 point) const
    {
        return findNearestPoints(point, 1)[0];
    };
    real distanceToPoint(uint index, real3 point) const
    {
        const real3 dist = coordinates[index] - point;
        return std::hypot(dist.x, dist.y, dist.z);
    };
    void printTree()
    {
        applyToTree(root.get(), [&](uint index, uint depth) {
            VF_LOG_INFO("index {} coords {},{},{} depth {}", index, coordinates[index].x, coordinates[index].y,
                        coordinates[index].z, depth);
        });
    }

private:
    void fillTree()
    {
        this->root = insertNodes(coordinates);
    }
    std::unique_ptr<Node> root;
    const std::vector<real3> coordinates;
    const std::vector<T> values;
    void nearestNeighborHelper(const Node* node, const real3& point, uint depth, uint& best, real& bestDist) const;
    const std::unique_ptr<InterpolationStrategy> interpolationStrategy;
};

template <typename T>
class GbPointCloud3D<T>::InterpolationStrategy
{
public:
    virtual T interpolateToPoint(real3 point) const = 0;
    void setData(GbPointCloud3D<T>* data)
    {
        this->data = data;
    };
    virtual ~InterpolationStrategy() {};
    GbPointCloud3D<T>* data;
};

template <typename T>
class GbPointCloud3D<T>::InverseDistanceWeighing : public GbPointCloud3D<T>::InterpolationStrategy
{
    using InterpolationStrategy = typename GbPointCloud3D<T>::InterpolationStrategy;

public:
    static std::unique_ptr<InterpolationStrategy> make(uint nPoints = 10)
    {
        return std::unique_ptr<InterpolationStrategy>(new InverseDistanceWeighing(nPoints));
    }
    InverseDistanceWeighing(uint nPoints = 10) : nPoints(nPoints) {};
    T interpolateToPoint(real3 point) const override
    {
        const auto nearestNeighbors = this->data->findNearestPoints(point, nPoints);
        real sumWeights {};
        T sumValues {};
        for (auto index : nearestNeighbors) {
            const real distance = this->data->distanceToPoint(index, point);
            if(distance < vf::basics::constant::cSmallSingle)
                return this->data->getValue(index);
            const real weight = std::pow(distance, -2);
            sumValues += weight * this->data->getValue(index);
            sumWeights += weight;
        }
        return sumValues / sumWeights;
    };

private:
    const uint nPoints;
};

template <typename T>
class GbPointCloud3D<T>::NearestNeighbor : public GbPointCloud3D<T>::InterpolationStrategy
{
    using InterpolationStrategy = typename GbPointCloud3D<T>::InterpolationStrategy;

public:
    static std::unique_ptr<InterpolationStrategy> make()
    {
        return std::unique_ptr<InterpolationStrategy>(new NearestNeighbor());
    }
    NearestNeighbor() = default;
    T interpolateToPoint(real3 point) const override
    {
        return this->data->getValue(this->data->findNearestPoint(point));
    };
};

#endif
