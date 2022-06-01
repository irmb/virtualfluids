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
//! \file WbWriterBOBJ.h
//! \ingroup writer
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifdef CAB_ZLIB
#ifndef WBWRITERBOBJ_H
#define WBWRITERBOBJ_H

#include <basics/writer/WbWriter.h>
#include <string>

class WbWriterBOBJ : public WbWriter
{
public:
    OBCREATOR_EXT(WbWriterBOBJ)

    static WbWriterBOBJ *getInstance()
    {
        static WbWriterBOBJ instance;
        return &instance;
    }

private:
    WbWriterBOBJ() : WbWriter()
    {
        if (sizeof(unsigned char) != 1)
            throw UbException(UB_EXARGS, "error char  type mismatch");
        if (sizeof(int) != 4)
            throw UbException(UB_EXARGS, "error int   type mismatch");
        if (sizeof(float) != 4)
            throw UbException(UB_EXARGS, "error float type mismatch");
    }
    WbWriterBOBJ(const WbWriterBOBJ &);                  // no copy allowed
    const WbWriterBOBJ &operator=(const WbWriterBOBJ &); // no copy allowed

    static std::string pvdEndTag;

public:
    std::string getFileExtension() { return "BOBJ.gz"; }

    std::string writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                               std::vector<UbTupleInt3> &triangles);
};

UB_AUTO_RUN_NAMED(
    ObFactory<WbWriter>::getInstance()->addObCreator(ObSingletonCreatorImpl<WbWriterBOBJ, WbWriter>::getInstance()),
    CAB_WbWriterVtkXmlASCII);

#endif // WBWRITERBOBJ_H

#endif // CAB_ZLIB
