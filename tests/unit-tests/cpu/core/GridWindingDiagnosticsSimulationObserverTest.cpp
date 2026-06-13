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
//! \addtogroup cpu_core_tests core
//! \ingroup cpu_core_tests
//! \{
//! \author Hussein Alihussein
//! \brief Validates missing-link diagnostics write a winding value for an open surface.
//!
//! Flow:
//! 1) Construct a Grid3D with 2 nodes in X, 2 in Y, and 2 in Z, and an open tetrahedron
//!    with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1). The triangles are
//!    (0,1,3), (0,3,2), (1,2,3), leaving the (0,1,2) face open.
//! 2) Build one missing link against the open surface: use the mesh triangles directly,
//!    then processBoundaryLink() traces a single link segment from origin
//!    (0.2, 0.2, -0.2) to neighbour (0.2, 0.2, 0.2) using dx=0.4 and (cx,cy,cz)=(0,0,1).
//!    The surface is open, so there is no hit;
//!    handleMissing() returns false and the link is recorded as missing.
//! 3) Finalize the collection: finalizeMissingLinkCollection() uses the fast-winding
//!    evaluator (when VF_HAS_FAST_WINDING is available) to compute and attach winding
//!    values for the missing-link endpoints, and also computes endpoint-kind flags.
//! 4) Write the collection to VTK and read back the winding scalar from the file.
//! 5) Parse the appended raw VTK data to extract the winding scalar and assert it is
//!    finite and within the expected [0, 1.1] tolerance band.
//!
//! Additional intersection cases covered below using the same tetrahedron mesh:
//! - Segment ending exactly on a triangle (endpoint hit at neighbour).
//! - Segment crossing a triangle in the middle (interior hit).
//! - Segment starting exactly on a triangle (endpoint hit at origin).
//! - Segment entirely outside all triangles (no hit, no-intersection link).
//=======================================================================================

#include <gtest/gtest.h>

#include <basics/PointerDefinitions.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/geometry3d/GbSystem3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>
#include <basics/geometry3d/winding/GridWindingWriting.h>
#include <cpu/core/SimulationObservers/GridWindingDiagnosticsSimulationObserver.h>
#include <parallel/NullCommunicator.h>

#include "Grid3D.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

namespace
{
SPtr<GbTriFaceMesh3D> makeOpenTetrahedron()
{
    auto mesh  = std::make_shared<GbTriFaceMesh3D>();
    auto *nodes = mesh->getNodes();
    auto *tris  = mesh->getTriangles();

    nodes->emplace_back(0.0f, 0.0f, 0.0f);
    nodes->emplace_back(1.0f, 0.0f, 0.0f);
    nodes->emplace_back(0.0f, 1.0f, 0.0f);
    nodes->emplace_back(0.0f, 0.0f, 1.0f);

    tris->emplace_back(0, 1, 3);
    tris->emplace_back(0, 3, 2);
    tris->emplace_back(1, 2, 3);

    mesh->calculateValues();
    return mesh;
}

vf::grid_winding::MissingLinkCollection buildMissingLinksForOpenSurface(
    const SPtr<GbTriFaceMesh3D> &surface,
    const SPtr<vf::parallel::Communicator> &comm)
{
    using vf::grid_winding::Vec3;
    vf::grid_winding::SubgridDistanceStats stats;
    std::vector<GbTriFaceMesh3D*> surfaces{ surface.get() };

    Vec3 origin(0.2, 0.2, -0.2);
    const double dx  = 0.4;
    const int    cx  = 0;
    const int    cy  = 0;
    const int    cz  = 1;
    const int    dir = 1;

    auto promoteSolid = []() { return false; };
    auto existingQ = []() { return -1.0; };
    auto setQ = [](double, int) {};
    auto markSurface = [](double, int) {};
    auto getWindingNumber = []() -> std::pair<bool, double> { return { false, 0.0 }; };
    auto handleMissing = [](int, const Vec3 &, const Vec3 &,
                            vf::grid_winding::SubgridDistanceStats::MissingLink &) { return false; };

    vf::grid_winding::processBoundaryLink(surfaces, origin, dx, cx, cy, cz, dir,
                                          promoteSolid, existingQ, setQ, markSurface,
                                          getWindingNumber, handleMissing, &stats);

    return vf::grid_winding::collectMissingLinks(stats, comm);
}

struct LinkTraceResult
{
    bool   hit{ false };
    bool   missing{ false };
    bool   qSet{ false };
    double q{ std::numeric_limits<double>::quiet_NaN() };
};

LinkTraceResult traceBoundaryLink(const SPtr<GbTriFaceMesh3D> &surface,
                                  const vf::grid_winding::Vec3 &origin,
                                  double dx, int cx, int cy, int cz)
{
    using vf::grid_winding::Vec3;
    vf::grid_winding::SubgridDistanceStats stats;
    LinkTraceResult result;
    std::vector<GbTriFaceMesh3D*> surfaces{ surface.get() };

    auto promoteSolid = []() { return false; };
    auto existingQ = []() { return -1.0; };
    auto setQ = [&](double q, int) {
        result.qSet = true;
        result.q = q;
    };
    auto markSurface = [](double, int) {};
    auto getWindingNumber = []() -> std::pair<bool, double> { return { false, 0.0 }; };
    auto handleMissing = [&](int, const Vec3 &, const Vec3 &,
                             vf::grid_winding::SubgridDistanceStats::MissingLink &) {
        result.missing = true;
        return false;
    };

    result.hit = vf::grid_winding::processBoundaryLink(surfaces, origin, dx, cx, cy, cz, 1,
                                                       promoteSolid, existingQ, setQ, markSurface,
                                                       getWindingNumber, handleMissing, &stats);
    if (!stats.missingLinks.empty())
        result.missing = true;

    return result;
}

std::filesystem::path makeTempOutputPath()
{
    const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::filesystem::path base = std::filesystem::temp_directory_path() / "vf_grid_winding_test";
    std::filesystem::path path = base / std::to_string(static_cast<long long>(now));
    std::filesystem::create_directories(path);
    return path;
}

bool readWindingFromFile(const std::filesystem::path &file, float &out)
{
    std::ifstream in(file, std::ios::binary);
    if (!in)
        return false;

    std::vector<char> data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (data.empty())
        return false;

    const std::string marker = "<AppendedData encoding=\"raw\">";
    auto markerIt = std::search(data.begin(), data.end(), marker.begin(), marker.end());
    if (markerIt == data.end())
        return false;

    auto underscoreIt = std::find(markerIt, data.end(), '_');
    if (underscoreIt == data.end())
        return false;

    std::size_t offset = static_cast<std::size_t>(std::distance(data.begin(), underscoreIt) + 1);

    auto readInt32 = [&](int32_t &value) -> bool {
        if (offset + sizeof(int32_t) > data.size())
            return false;
        std::memcpy(&value, &data[offset], sizeof(int32_t));
        offset += sizeof(int32_t);
        return true;
    };

    auto readFloat = [&](float &value) -> bool {
        if (offset + sizeof(float) > data.size())
            return false;
        std::memcpy(&value, &data[offset], sizeof(float));
        offset += sizeof(float);
        return true;
    };

    auto skipBytes = [&](std::size_t count) -> bool {
        if (offset + count > data.size())
            return false;
        offset += count;
        return true;
    };

    int32_t bytes = 0;
    if (!readInt32(bytes) || bytes < 0 || !skipBytes(static_cast<std::size_t>(bytes)))
        return false;
    if (!readInt32(bytes) || bytes < 0 || !skipBytes(static_cast<std::size_t>(bytes)))
        return false;
    if (!readInt32(bytes) || bytes < 0 || !skipBytes(static_cast<std::size_t>(bytes)))
        return false;
    if (!readInt32(bytes) || bytes < 0 || !skipBytes(static_cast<std::size_t>(bytes)))
        return false;

    if (!readInt32(bytes) || bytes != static_cast<int32_t>(sizeof(float)))
        return false;
    float dir = 0.0f;
    if (!readFloat(dir))
        return false;
    (void)dir;

    if (!readInt32(bytes) || bytes != static_cast<int32_t>(sizeof(float)))
        return false;
    float endpoint = 0.0f;
    if (!readFloat(endpoint))
        return false;
    (void)endpoint;

    if (!readInt32(bytes) || bytes != static_cast<int32_t>(sizeof(float)))
        return false;
    return readFloat(out);
}

void writeVisualizationGeometry(const SPtr<GbTriFaceMesh3D> &surface,
                                const std::filesystem::path &outputPath,
                                const std::filesystem::path &missingLinkFile)
{
    std::filesystem::path vizPath = outputPath / "viz";
    std::filesystem::create_directories(vizPath);

    const std::string surfaceFile = gb_system_3d::writeGeoObject(
        surface,
        (vizPath / "tetrahedron_surface").string(),
        WbWriterVtkXmlBinary::getInstance());
    std::cerr << "Tetrahedron surface VTK file: " << surfaceFile << std::endl;

    struct SegmentCase
    {
        const char *name;
        vf::grid_winding::Vec3 origin;
        vf::grid_winding::Vec3 neighbour;
        float expectedHit;
        float expectedMissing;
    };

    std::vector<SegmentCase> cases;
    cases.push_back({ "missing_link",
                      vf::grid_winding::Vec3(0.2, 0.2, -0.2),
                      vf::grid_winding::Vec3(0.2, 0.2, 0.2),
                      0.0f,
                      1.0f });
    cases.push_back({ "hit_neighbour_endpoint",
                      vf::grid_winding::Vec3(-0.2, 0.1, 0.1),
                      vf::grid_winding::Vec3(0.0, 0.1, 0.1),
                      1.0f,
                      0.0f });
    cases.push_back({ "hit_middle",
                      vf::grid_winding::Vec3(-0.2, 0.3, 0.1),
                      vf::grid_winding::Vec3(0.2, 0.3, 0.1),
                      1.0f,
                      0.0f });
    cases.push_back({ "hit_origin_endpoint",
                      vf::grid_winding::Vec3(0.0, 0.6, 0.1),
                      vf::grid_winding::Vec3(0.2, 0.6, 0.1),
                      1.0f,
                      0.0f });
    cases.push_back({ "outside_no_hit",
                      vf::grid_winding::Vec3(2.0, 2.0, 2.0),
                      vf::grid_winding::Vec3(2.4, 2.0, 2.0),
                      0.0f,
                      1.0f });

    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2> lines;
    nodes.reserve(cases.size() * 2);
    lines.reserve(cases.size());

    std::vector<std::string> lineDataNames{ "case_id", "expected_hit", "expected_missing" };
    std::vector<std::vector<float>> lineData(3);
    lineData[0].reserve(cases.size());
    lineData[1].reserve(cases.size());
    lineData[2].reserve(cases.size());

    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto &entry = cases[i];
        const int   nodeIndex = static_cast<int>(nodes.size());
        nodes.emplace_back(static_cast<float>(entry.origin.X1()),
                           static_cast<float>(entry.origin.X2()),
                           static_cast<float>(entry.origin.X3()));
        nodes.emplace_back(static_cast<float>(entry.neighbour.X1()),
                           static_cast<float>(entry.neighbour.X2()),
                           static_cast<float>(entry.neighbour.X3()));
        lines.emplace_back(nodeIndex, nodeIndex + 1);
        lineData[0].push_back(static_cast<float>(i));
        lineData[1].push_back(entry.expectedHit);
        lineData[2].push_back(entry.expectedMissing);
        std::cerr << "Segment case " << i << ": " << entry.name << std::endl;
    }

    const std::string segmentsFile = WbWriterVtkXmlBinary::getInstance()->writeLinesWithLineData(
        (vizPath / "boundary_link_cases").string(),
        nodes,
        lines,
        lineDataNames,
        lineData);
    std::cerr << "Boundary link cases VTK file: " << segmentsFile << std::endl;

    const std::filesystem::path manifest = vizPath / "paths.txt";
    std::ofstream out(manifest);
    if (out) {
        out << "missing_links=" << missingLinkFile << "\n";
        out << "tetrahedron_surface=" << surfaceFile << "\n";
        out << "boundary_link_cases=" << segmentsFile << "\n";
    }
    std::cerr << "Visualization manifest: " << manifest << std::endl;
}
} // namespace

TEST(GridWindingDiagnosticsSimulationObserverTest, WritesMissingLinkWinding)
{
#if !defined(VF_HAS_FAST_WINDING)
    GTEST_SKIP() << "Fast winding support not available.";
#endif

    auto comm = vf::parallel::NullCommunicator::getInstance();
    auto grid = std::make_shared<Grid3D>(comm, 2, 2, 2, 1, 1, 1);
    auto surface = makeOpenTetrahedron();

    const auto outputPath = makeTempOutputPath();
    auto collection = buildMissingLinksForOpenSurface(surface, comm);
    ASSERT_FALSE(collection.empty());
    vf::grid_winding::cpu::finalizeMissingLinkCollection(collection, surface, grid);
    ASSERT_FALSE(collection.empty());

    std::filesystem::create_directories(outputPath / "bc");
    vf::grid_winding::writeMissingLinks(collection, (outputPath / "bc" / "missing_boundary_links").string());

    const std::filesystem::path file = outputPath / "bc" / "missing_boundary_links.bin.vtu";
    std::cerr << "Missing-link VTK file: " << file << std::endl;
    writeVisualizationGeometry(surface, outputPath, file);
    ASSERT_TRUE(std::filesystem::exists(file));

    float winding = std::numeric_limits<float>::quiet_NaN();
    ASSERT_TRUE(readWindingFromFile(file, winding));
    ASSERT_TRUE(std::isfinite(winding));
    EXPECT_TRUE(ub_math::greaterEqual(winding, 0.0f));
    EXPECT_TRUE(ub_math::lessEqual(winding, 1.1f));
}

TEST(GridWindingDiagnosticsSimulationObserverTest, BoundaryLinkHitsAtNeighbourEndpoint)
{
    auto surface = makeOpenTetrahedron();
    const vf::grid_winding::Vec3 origin(-0.2, 0.1, 0.1);

    const auto result = traceBoundaryLink(surface, origin, 0.2, 1, 0, 0);
    EXPECT_TRUE(result.hit);
    EXPECT_FALSE(result.missing);
    ASSERT_TRUE(result.qSet);
    EXPECT_NEAR(result.q, 1.0, 1.0e-6);
}

TEST(GridWindingDiagnosticsSimulationObserverTest, BoundaryLinkHitsInMiddle)
{
    auto surface = makeOpenTetrahedron();
    const vf::grid_winding::Vec3 origin(-0.2, 0.3, 0.1);

    const auto result = traceBoundaryLink(surface, origin, 0.4, 1, 0, 0);
    EXPECT_TRUE(result.hit);
    EXPECT_FALSE(result.missing);
    ASSERT_TRUE(result.qSet);
    EXPECT_NEAR(result.q, 0.5, 1.0e-6);
}

TEST(GridWindingDiagnosticsSimulationObserverTest, BoundaryLinkHitsAtOriginEndpoint)
{
    auto surface = makeOpenTetrahedron();
    const vf::grid_winding::Vec3 origin(0.0, 0.6, 0.1);

    const auto result = traceBoundaryLink(surface, origin, 0.2, 1, 0, 0);
    EXPECT_TRUE(result.hit);
    EXPECT_FALSE(result.missing);
    ASSERT_TRUE(result.qSet);
    EXPECT_NEAR(result.q, 0.0, 1.0e-6);
}

TEST(GridWindingDiagnosticsSimulationObserverTest, BoundaryLinkOutsideSegmentNoHit)
{
    auto surface = makeOpenTetrahedron();
    const vf::grid_winding::Vec3 origin(2.0, 2.0, 2.0);

    const auto result = traceBoundaryLink(surface, origin, 0.4, 1, 0, 0);
    EXPECT_FALSE(result.hit);
    EXPECT_TRUE(result.missing);
    EXPECT_FALSE(result.qSet);
}

//! \}
