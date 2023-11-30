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
//! \file STLWriter.cpp
//! \ingroup io
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include "STLWriter.h"

#include <fstream>
#include <sstream>

#include "geometries/Vertex/Vertex.h"
#include "geometries/Triangle/Triangle.h"

void STLWriter::writeSTL(std::vector<Triangle> &vec, const std::string &name, bool writeBinary)
{
    const int size = (int)vec.size();
    VF_LOG_INFO("Write {} Triangles to STL : {}" , size, name);

    std::ios_base::openmode mode = std::ios::out;
    if (writeBinary)
        mode = std::ios::out | std::ios::binary;

    std::ofstream ofstream(name, mode);

    if (!ofstream.is_open()) {
        VF_LOG_WARNING("Output file not open - exit function");
        return;
    }

    if (writeBinary)
        writeBinarySTL(ofstream, vec);
    else
        writeAsciiSTL(ofstream, vec);

    ofstream.close();
    VF_LOG_INFO("Output file closed");
}


void STLWriter::writeAsciiSTL(std::ofstream &ofstream, std::vector<Triangle> &vec)
{
    ofstream << "solid ascii\n";
    for (size_t i = 0; i < vec.size(); i++) 
    {
        Triangle t = vec[i];

        ofstream << "facet normal ";
        t.normal.printFormatted(ofstream);
        ofstream << "\n";

        ofstream << "outer loop\n";

        ofstream << "vertex ";
        t.v1.printFormatted(ofstream);
        ofstream << "\n";
        ofstream << "vertex ";
        t.v2.printFormatted(ofstream);
        ofstream << "\n";
        ofstream << "vertex ";
        t.v3.printFormatted(ofstream);
        ofstream << "\n";

        ofstream << "endloop\n";
        ofstream << "endfacet\n";
    }
    ofstream << "endsolid\n";
}

void STLWriter::writeBinarySTL(std::ofstream &ofstream, std::vector<Triangle> &vec)
{
    char header_info[80] = "GridGeneration-File iRMB";
    unsigned long nTriLong = (unsigned long)vec.size();
    ofstream.write(header_info, sizeof(header_info));
    ofstream.write((char*)&nTriLong, 4);
    char attribute[2] = "0";
    for (size_t i = 0; i < vec.size(); i++)
    {
        Triangle t = vec[i];

        t.normal.print(ofstream);

        t.v1.print(ofstream);
        t.v2.print(ofstream);
        t.v3.print(ofstream);

        ofstream.write(attribute, 2);
    }
}
