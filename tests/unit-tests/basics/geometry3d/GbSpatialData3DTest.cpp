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
//! \addtogroup utilities_test utilities
//! \ingroup basics_test
//! \{
//! \author Henry Korb
//=======================================================================================
#include "DataTypes.h"
#include "constants/NumericConstants.h"
#include "gmock/gmock.h"

#include "gtest/gtest.h"
#include <basics/geometry3d/GbSpatialData3D.h>
#include <cmath>
#include <memory>

using namespace testing;
using namespace vf::basics::constant;

const real c1o2p25 = c2o3 * c2o3;
template<typename T> internal::FloatingEqMatcher<T> Equal(T a);
template <> internal::FloatingEqMatcher<float> Equal(float a){return FloatEq(a);};
template <> internal::FloatingEqMatcher<double> Equal(double a){return DoubleEq(a);};


TEST(StructuredMeshTest, Interpolation3DReal)
{
    const real3 spacing { 1.0, 1.0, 1.0 };
    const real3 origin {};
    const std::array<uint, 3> nPoints { 2, 2, 2 };
    const std::vector<real> data { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
    GbStructuredMesh3D<real> mesh(spacing, origin, nPoints, data,
                                  GbStructuredMesh3D<real>::InterpolationStrategy::Trilinear);
    EXPECT_THAT(mesh.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(0.5));
    EXPECT_THAT(mesh.interpolateToPoint({ 0, 0., 1.0 }), Eq(1.0));
}

TEST(StructuredMeshTest, Interpolation3DReal3)
{
    const real3 spacing { 1.0, 1.0, 1.0 };
    const real3 origin {};
    const std::array<uint, 3> nPoints { 2, 2, 2 };
    const real3 zero {};
    const real3 one { 1.0, 1.0, 1.0 };
    const std::vector<real3> data { zero, zero, zero, zero, one, one, one, one };
    GbStructuredMesh3D<real3> mesh(spacing, origin, nPoints, data,
                                   GbStructuredMesh3D<real3>::InterpolationStrategy::Trilinear);
    EXPECT_THAT(mesh.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(real3 { 0.5, 0.5, 0.5 }));
    EXPECT_THAT(mesh.interpolateToPoint({ 0, 0., 1.0 }), Eq(one));
}

TEST(StructuredMeshTest, Interpolation2DReal)
{
    const real3 spacing { 1.0, 1.0, 0.0 };
    const real3 origin {};
    const std::array<uint, 3> nPoints { 2, 2, 1 };
    const std::vector<real> data { 0.0, 1.0, 0.0, 1.0 };
    GbStructuredMesh3D<real> mesh(spacing, origin, nPoints, data,
                                  GbStructuredMesh3D<real>::InterpolationStrategy::Trilinear);
    EXPECT_THAT(mesh.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(0.5));
    EXPECT_THAT(mesh.interpolateToPoint({ 0, 0., 1.0 }), Eq(0.0));
    EXPECT_THAT(mesh.interpolateToPoint({ 1.0, 1.0, -1.5 }), Eq(1.0));
}

TEST(PointCloudTest, InterpolationRealIDW3D)
{
    std::vector<real3> points {
        {},
    };
    std::vector<real> data { 1.0 };
    GbPointCloud3D<real> pointCloud(points, data, std::make_unique<GbPointCloud3D<real>::InverseDistanceWeighing>());
    EXPECT_THAT(pointCloud.interpolateToPoint({}), Eq(1.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(1.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, -0.5, -0.5 }), Eq(1.0));
}

TEST(PointCloudTest, InterpolationRealIDW3D2)
{
    std::vector<real3> points { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };
    std::vector<real> data { 0.0, 1.0 };
    GbPointCloud3D<real> pointCloud(points, data, std::make_unique<GbPointCloud3D<real>::InverseDistanceWeighing>());
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.0, 0.0, 0.0 }), Eq(0.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(0.5));
    const real weights = c4o1 + c1o2p25;
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, 0.0, 0.0 }), Equal(c1o2p25 / weights));
}

TEST(PointCloudTest, InterpolationRealIDW3D3)
{
    std::vector<real3> points { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 } };
    std::vector<real> data { 0.0, 1.0, 100.0 };
    GbPointCloud3D<real> pointCloud(points, data, std::make_unique<GbPointCloud3D<real>::InverseDistanceWeighing>(2));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.0, 0.0, 0.0 }), Eq(0.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(0.5));
    const real weights = c4o1 + c1o2p25;
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, 0.0, 0.0 }), Equal(c1o2p25 / weights));
}

TEST(PointCloudTest, InterpolationReal3IDW3D)
{
    std::vector<real3> points {
        {},
    };
    const real3 one { 1.0, 1.0, 1.0 };
    std::vector<real3> data { one };
    GbPointCloud3D<real3> pointCloud(points, data, std::make_unique<GbPointCloud3D<real3>::InverseDistanceWeighing>());
    EXPECT_THAT(pointCloud.interpolateToPoint({}), Eq(one));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(one));
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, -0.5, -0.5 }), Eq(one));
}

TEST(PointCloudTest, InterpolationReal3IDW3D2)
{
    std::vector<real3> points { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };
    const real3 zero {};
    const real3 one { 1.0, 1.0, 1.0 };
    std::vector<real3> data { zero, one };
    GbPointCloud3D<real3> pointCloud(points, data, GbPointCloud3D<real3>::InverseDistanceWeighing::make());
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.0, 0.0, 0.0 }), Eq(zero));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(one * 0.5));
    const real weights = c4o1 + c1o2p25;
    const real3 result = pointCloud.interpolateToPoint({ -0.5, 0.0, 0.0 });
    EXPECT_THAT(result.x, Equal(c1o2p25 / weights));
    EXPECT_THAT(result.y, Equal(c1o2p25 / weights));
    EXPECT_THAT(result.z, Equal(c1o2p25 / weights));
}

TEST(PointCloudTest, InterpolationNearestNeighborReal3D)
{
    std::vector<real3> points {
        {},
    };
    std::vector<real> data { 1.0 };
    GbPointCloud3D<real> pointCloud(points, data, GbPointCloud3D<real>::NearestNeighbor::make());
    EXPECT_THAT(pointCloud.interpolateToPoint({}), Eq(1.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(1.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, -0.5, -0.5 }), Eq(1.0));
}

TEST(PointCloudTest, InterpolationNearestNeighborReal3D2)
{
    std::vector<real3> points { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };
    std::vector<real> data { 0.0, 1.0 };
    GbPointCloud3D<real> pointCloud(points, data, GbPointCloud3D<real>::NearestNeighbor::make());
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.0, 0.0, 0.0 }), Eq(0.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(1.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.49999, 0.5, 0.5 }), Eq(0.0));
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, 0.0, 0.0 }), Eq(0.0));
}

TEST(PointCloudTest, InterpolationNearestNeighborReal3D3)
{
    std::vector<real3> points {
        {},
    };
    const real3 one { 1.0, 1.0, 1.0 };
    std::vector<real3> data { one };
    GbPointCloud3D<real3> pointCloud(points, data, GbPointCloud3D<real3>::NearestNeighbor::make());
    EXPECT_THAT(pointCloud.interpolateToPoint({}), Eq(one));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(one));
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, -0.5, -0.5 }), Eq(one));
}

TEST(PointCloudTest, InterpolationNearestNeighborReal33D2)
{
    std::vector<real3> points { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };
    const real3 zero {};
    const real3 one { 1.0, 1.0, 1.0 };
    std::vector<real3> data { zero, one };
    GbPointCloud3D<real3> pointCloud(points, data, GbPointCloud3D<real3>::NearestNeighbor::make());
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.0, 0.0, 0.0 }), Eq(zero));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.5, 0.5, 0.5 }), Eq(one));
    EXPECT_THAT(pointCloud.interpolateToPoint({ 0.49999, 0.5, 0.5 }), Eq(zero));
    EXPECT_THAT(pointCloud.interpolateToPoint({ -0.5, 0.0, 0.0 }), Eq(zero));
}