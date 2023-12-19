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
//! \addtogroup gpu_io io
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef SimulationFileWriter_H
#define SimulationFileWriter_H

#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <vector>
#include <array>

#include "gpu/GridGenerator/global.h"

class UnstructuredGridBuilder;
class GridBuilder;
class Grid;
namespace gg
{
class BoundaryCondition;
}

enum class FILEFORMAT
{
    BINARY, ASCII
};

class SimulationFileWriter
{
public:
    static void write(const std::string& folder, SPtr<GridBuilder> builder, FILEFORMAT format);

private:
    static void write(SPtr<GridBuilder> builder, FILEFORMAT format);
    static void openFiles(SPtr<GridBuilder> builder);

    static void writeNumberNodes(SPtr<GridBuilder> builder, uint level);
    static void writeLBMvsSI(SPtr<GridBuilder> builder, uint level);

    static void writeLevel(uint numberOfLevels);
    static void writeLevelSize(uint numberOfNodes, FILEFORMAT format);
    static void writeCoordFiles(SPtr<GridBuilder> builder, uint level, FILEFORMAT format);
    static void writeCoordsNeighborsGeo(SPtr<GridBuilder> builder, int index, uint level, FILEFORMAT format);

    static void writeLevelSizeGridInterface(uint sizeCF, uint sizeFC);
    static void writeGridInterfaceToFile(SPtr<GridBuilder> builder, uint level);
    static void writeGridInterfaceToFile(uint numberOfNodes, std::ofstream& coarseFile, uint* coarse, std::ofstream& fineFile, uint* fine);
    static void writeGridInterfaceOffsetToFile(uint numberOfNodes, std::ofstream& offsetFile, real* cf_offset_X, real* cf_offset_Y, real* cf_offset_Z);

    static void writeBoundaryQsFile(SPtr<GridBuilder> builder);
    static std::vector<std::vector<std::vector<real> > > createBCVectors(SPtr<Grid> grid);
    static void addShortQsToVector(int index, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid);
    static void addQsToVector(int index, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid);
    static void fillRBForNode(int index, int direction, int directionSign, int rb, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid);
    static void writeBoundary(std::vector<real> boundary, int rb);
    static void writeBoundaryShort(std::vector<real> boundary, int rb);
    static void writeBoundaryShort(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundaryCondition, uint side);

    static void writeCommunicationFiles(SPtr<GridBuilder> builder);

    static void closeFiles();


    static std::ofstream xCoordFile;
    static std::ofstream yCoordFile;
    static std::ofstream zCoordFile;
    static std::ofstream xNeighborFile;
    static std::ofstream yNeighborFile;
    static std::ofstream zNeighborFile;
    static std::ofstream wsbNeighborFile;
    static std::ofstream geoVecFile;

    static std::ofstream scaleCF_coarse_File;
    static std::ofstream scaleCF_fine_File;
    static std::ofstream scaleFC_coarse_File;
    static std::ofstream scaleFC_fine_File;

    static std::ofstream offsetVecCF_File;
    static std::ofstream offsetVecFC_File;

    static std::vector<SPtr<std::ofstream> > qStreams;
    static std::vector<SPtr<std::ofstream> > valueStreams;

    static std::vector<std::ofstream> qFiles;
    static std::vector<std::ofstream> valueFiles;

    static std::array<std::ofstream, 6> sendFiles;
    static std::array<std::ofstream, 6> receiveFiles;

    static std::ofstream numberNodes_File;
    static std::ofstream LBMvsSI_File;

    static std::string folder;
};


#endif

//! \}
