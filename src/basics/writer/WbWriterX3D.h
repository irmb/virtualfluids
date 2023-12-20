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
//! \addtogroup writer
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef WBWRITERX3D_H
#define WBWRITERX3D_H

#include <string>

#include <basics/writer/WbWriter.h>

class WbWriterX3D : public WbWriter
{
public:
    static WbWriterX3D *getInstance()
    {
        static WbWriterX3D instance;
        return &instance;
    }

    WbWriterX3D(const WbWriterX3D &) = delete;
    const WbWriterX3D &operator=(const WbWriterX3D &) = delete;

private:
    WbWriterX3D() : WbWriter()
    {
        if constexpr (sizeof(unsigned char) != 1)
            throw UbException(UB_EXARGS, "error char  type mismatch");
        if constexpr (sizeof(int) != 4)
            throw UbException(UB_EXARGS, "error int   type mismatch");
        if constexpr (sizeof(float) != 4)
            throw UbException(UB_EXARGS, "error float type mismatch");
    }

    static std::string pvdEndTag;

public:
    std::string getFileExtension() override { return "ascii.X3D"; }

    std::string writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                               std::vector<UbTupleInt3> &triangles) override;
};

#endif // WBWRITERX3D_H

//! \}
