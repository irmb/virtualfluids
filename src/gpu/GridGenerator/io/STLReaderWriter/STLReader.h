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
#ifndef STLReader_H
#define STLReader_H


#include <vector>
#include <string>

#include "global.h"

struct Triangle;
struct Vertex;
class BoundingBox;

class STLReader
{
public:

    enum FileType { ascii, binary };

    static std::vector<Triangle> readSTL(const std::string& name);
    static std::vector<Triangle> readSTL(const std::string& name, FileType fileType, const std::vector<uint> ignorePatches = std::vector<uint>() );
    static std::vector<Triangle> readSTL(const BoundingBox &box, const std::string& name);

    static std::vector<Triangle> readBinarySTL(const std::string& name);
    static std::vector<Triangle> readASCIISTL(const std::string& name);
    static std::vector<Triangle> readASCIISTLWithPatches(const std::string& name, const std::vector<uint> ignorePatches);
    static std::vector<Triangle> readBinarySTL(const BoundingBox &box, const std::string& name);
    static std::vector<Triangle> readASCIISTL(const BoundingBox &box, const std::string& name);

private:
    STLReader(){};
    ~STLReader(){};

    static int countLinesInFile(std::string name);
    static Vertex parseLineToCoordinates(std::ifstream& file, std::string format);
    static Vertex parseLineToCoordinates(const std::string& file, const std::string format);
    static Vertex getVertexFromChar(const char* facet);
};

#endif

//! \}
