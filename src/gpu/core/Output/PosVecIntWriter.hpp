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
#ifndef POSITION_VECTOR_INTEGER_WRITER_H
#define POSITION_VECTOR_INTEGER_WRITER_H

#include <basics/utilities/UbFileOutputASCII.h>

#include "Parameter/Parameter.h"

class PositionVectorIntegerWriter
{
public:
    static void writePositionInterface(Parameter* para, std::string Type)
    {
        UbFileOutputASCII out(para->getFName() + Type + ".dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        if (Type == "_InterfaceCFC") {
            for (int level = 0; level < para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->coarseToFine.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
                    out.writeInteger(para->getParH(level)->coarseToFine.coarseCellIndices[u]);
                }
                out.writeLine();
            }
        } else if (Type == "_InterfaceCFF") {
            for (int level = 0; level < para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->coarseToFine.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
                    out.writeInteger(para->getParH(level)->coarseToFine.fineCellIndices[u]);
                }
                out.writeLine();
            }
        } else if (Type == "_InterfaceFCC") {
            for (int level = 0; level < para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->fineToCoarse.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++) {
                    out.writeInteger(para->getParH(level)->fineToCoarse.coarseCellIndices[u]);
                }
                out.writeLine();
            }
        } else if (Type == "_InterfaceFCF") {
            for (int level = 0; level < para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->fineToCoarse.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++) {
                    out.writeInteger(para->getParH(level)->fineToCoarse.fineCellIndices[u]);
                }
                out.writeLine();
            }
        }
    }
};
#endif
//! \}
