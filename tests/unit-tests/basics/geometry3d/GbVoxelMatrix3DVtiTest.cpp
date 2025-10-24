
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
//! \author Konstantin Kutscher
//! \brief Unit tests for GbVoxelMatrix3D VTI file reading functionality. Generated with the assistance of GitHub Copilot ver. 1.387.0 using GPT-5 mini. Reviewed and adapted by Konstantin Kutscher.
//=======================================================================================

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdint>
#include <ctime>
#include <unistd.h>

#include <basics/geometry3d/GbVoxelMatrix3D.h>

static std::string makeTempFilename(const std::string &prefix)
{
    std::ostringstream ss;
    ss << "/tmp/" << prefix << "_" << getpid() << "_" << std::time(nullptr) << ".vti";
    return ss.str();
}

static void writeVtiAscii(const std::string &filename, int nx, int ny, int nz, const std::vector<float> &data)
{
    std::ofstream out(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(out.is_open());

    std::ostringstream xml;
    xml << "<?xml version=\"1.0\"?>\n";
    xml << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    xml << "<ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1)
        << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    xml << "  <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1) << "\">\n";
    xml << "    <PointData Scalars=\"scalars\">\n";
    xml << "      <DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\">\n";
    for (size_t i = 0; i < data.size(); ++i) {
        xml << data[i];
        if (i + 1 < data.size()) xml << " ";
    }
    xml << "\n      </DataArray>\n";
    xml << "    </PointData>\n";
    xml << "    <CellData/>\n";
    xml << "  </Piece>\n";
    xml << "</ImageData>\n";
    xml << "</VTKFile>\n";

    out << xml.str();
    out.close();
}

static void writeVtiAppended(const std::string &filename, int nx, int ny, int nz, const std::vector<float> &data)
{
    std::ofstream out(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(out.is_open());

    std::ostringstream xml;
    xml << "<?xml version=\"1.0\"?>\n";
    xml << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    xml << "<ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1)
        << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    xml << "  <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 " << (nz - 1) << "\">\n";
    xml << "    <PointData Scalars=\"scalars\">\n";
    // offset 0 -> first block directly after '_' (reader expects offset relative to underscore+1)
    xml << "      <DataArray type=\"Float32\" Name=\"scalars\" format=\"appended\" offset=\"0\"/>\n";
    xml << "    </PointData>\n";
    xml << "    <CellData/>\n";
    xml << "  </Piece>\n";
    xml << "</ImageData>\n";
    xml << "<AppendedData encoding=\"raw\">\n";

    out << xml.str();

    // VTK appended data starts with '_' char
    out.put('_');

    // block size (uint32) then raw bytes
    uint32_t rawSize = static_cast<uint32_t>(data.size() * sizeof(float));
    // write little-endian uint32
    uint8_t sizeBytes[4];
    sizeBytes[0] = static_cast<uint8_t>((rawSize >> 0) & 0xFF);
    sizeBytes[1] = static_cast<uint8_t>((rawSize >> 8) & 0xFF);
    sizeBytes[2] = static_cast<uint8_t>((rawSize >> 16) & 0xFF);
    sizeBytes[3] = static_cast<uint8_t>((rawSize >> 24) & 0xFF);
    out.write(reinterpret_cast<const char *>(sizeBytes), 4);

    // write float data in native representation (reader expects float bytes)
    for (size_t i = 0; i < data.size(); ++i) {
        float v = data[i];
        out.write(reinterpret_cast<const char *>(&v), sizeof(float));
    }

    out << "\n</AppendedData>\n";
    out << "</VTKFile>\n";
    out.close();
}

TEST(GbVoxelMatrix3DVtiTest, ReadMatrixFromVtiASCIIFile)
{
    const int nx = 2, ny = 2, nz = 1;
    // x fastest order: (0,0,0),(1,0,0),(0,1,0),(1,1,0)
    std::vector<float> data = {0.5f, 2.0f, 2.0f, 2.0f};

    std::string fname = makeTempFilename("test_ascii_vti");
    writeVtiAscii(fname, nx, ny, nz, data);

    GbVoxelMatrix3D sut(nx, ny, nz, GbVoxelMatrix3D::FLUID);
    sut.setThreshold(0.0, 1.0);

    ASSERT_NO_THROW(sut.readMatrixFromVtiASCIIFile(fname));

    EXPECT_FLOAT_EQ(sut(0, 0, 0), GbVoxelMatrix3D::SOLID);
    EXPECT_FLOAT_EQ(sut(1, 0, 0), GbVoxelMatrix3D::FLUID);
    EXPECT_FLOAT_EQ(sut(0, 1, 0), GbVoxelMatrix3D::FLUID);
    EXPECT_FLOAT_EQ(sut(1, 1, 0), GbVoxelMatrix3D::FLUID);

    std::remove(fname.c_str());
}

TEST(GbVoxelMatrix3DVtiTest, ReadMatrixFromVtiAppendedFile)
{
    const int nx = 2, ny = 2, nz = 1;
    std::vector<float> data = {0.5f, 2.0f, 2.0f, 2.0f};

    std::string fname = makeTempFilename("test_appended_vti");
    writeVtiAppended(fname, nx, ny, nz, data);

    GbVoxelMatrix3D sut(nx, ny, nz, GbVoxelMatrix3D::FLUID);
    sut.setThreshold(0.0, 1.0);

    ASSERT_NO_THROW(sut.readMatrixFromVtiAppendedFile(fname));

    EXPECT_FLOAT_EQ(sut(0, 0, 0), GbVoxelMatrix3D::SOLID);
    EXPECT_FLOAT_EQ(sut(1, 0, 0), GbVoxelMatrix3D::FLUID);
    EXPECT_FLOAT_EQ(sut(0, 1, 0), GbVoxelMatrix3D::FLUID);
    EXPECT_FLOAT_EQ(sut(1, 1, 0), GbVoxelMatrix3D::FLUID);

    std::remove(fname.c_str());
}