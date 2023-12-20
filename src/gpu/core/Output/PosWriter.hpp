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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#ifndef POSITION_WRITER_H
#define POSITION_WRITER_H

#include <basics/utilities/UbFileOutputASCII.h>

#include "Parameter/Parameter.h"

class PositionWriter
{
public:
    static void writePosition(Parameter* para, std::string Type)
    {
        UbFileOutputASCII out(para->getFName() + Type + ".dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        if (Type == "_geoSP") {
            for (int level = 0; level <= para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
                    out.writeInteger(para->getParH(level)->typeOfGridNode[index]);
                }
                out.writeLine();
            }
        } else if (Type == "_neighborX_SP") {
            for (int level = 0; level <= para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
                    out.writeInteger(para->getParH(level)->neighborX[index]);
                }
                out.writeLine();
            }
        } else if (Type == "_neighborY_SP") {
            for (int level = 0; level <= para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
                    out.writeInteger(para->getParH(level)->neighborY[index]);
                }
                out.writeLine();
            }
        } else if (Type == "_neighborZ_SP") {
            for (int level = 0; level <= para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
                    out.writeInteger(para->getParH(level)->neighborZ[index]);
                }
                out.writeLine();
            }
        }
    }
};
#endif

//! \}
