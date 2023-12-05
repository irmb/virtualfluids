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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef OFFSET_WRITER_H
#define OFFSET_WRITER_H

#include <basics/utilities/UbFileOutputASCII.h>

#include "Parameter/Parameter.h"

class OffsetWriter
{
public:
    static void writeOffset(Parameter* para, std::string Type)
    {
        UbFileOutputASCII out(para->getFName() + Type + ".dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        if (Type == "_OffsetCF") {
            for (int level = 0; level < para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->coarseToFine.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
                    out.writeDouble(para->getParH(level)->neighborCoarseToFine.x[u]);
                    out.writeDouble(para->getParH(level)->neighborCoarseToFine.y[u]);
                    out.writeDouble(para->getParH(level)->neighborCoarseToFine.z[u]);
                }
                out.writeLine();
            }
        } else if (Type == "_OffsetFC") {
            for (int level = 0; level < para->getMaxLevel(); level++) {
                out.writeInteger(para->getParH(level)->fineToCoarse.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++) {
                    out.writeDouble(para->getParH(level)->neighborFineToCoarse.x[u]);
                    out.writeDouble(para->getParH(level)->neighborFineToCoarse.y[u]);
                    out.writeDouble(para->getParH(level)->neighborFineToCoarse.z[u]);
                }
                out.writeLine();
            }
        }
    }
};

#endif
