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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file GridVTKWriter.h
//! \ingroup io
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef GridVTKWriter_h
#define GridVTKWriter_h

#include <string>

#include "global.h"

enum class WRITING_FORMAT { BINARY, ASCII };

class Grid;

class GridVTKWriter
{
public:
    static void writeSparseGridToVTK(SPtr<Grid> grid, const std::string& name, WRITING_FORMAT format = WRITING_FORMAT::ASCII);
    static void writeGridToVTKXML(SPtr<Grid> grid, const std::string& name);
    static void writeInterpolationCellsToVTKXML(SPtr<Grid> grid, SPtr<Grid> gridCoarse, const std::string& name);

private:
    GridVTKWriter() = default;
    ~GridVTKWriter() = default;

    static FILE *file;
    static WRITING_FORMAT format;

    static void initalVtkWriter(WRITING_FORMAT format, const std::string& name);

    static bool isBinaryWritingFormat();

    static void writeVtkFile(SPtr<Grid> grid);

    static void openFile(const std::string& name, const std::string& mode);
    static void closeFile();

    static void writeHeader();
    static void writePoints(SPtr<Grid> grid);
    static void writeCells(const unsigned int &size);
    static void writeTypeHeader(const unsigned int &size);
    static void writeTypes(SPtr<Grid> grid);

    static void end_line();
    static void force_big_endian(unsigned char *bytes);
    static void write_int(int val);
    static void write_float(float val);
};


#endif
