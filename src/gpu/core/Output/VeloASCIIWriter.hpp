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
#ifndef VELO_ASCII_WRITER_H
#define VELO_ASCII_WRITER_H

#include <numeric>

#include <basics/StringUtilities/StringUtil.h>
#include <basics/utilities/UbFileOutputASCII.h>

#include "Parameter/Parameter.h"

class VeloASCIIWriter
{
public:
    static void writeVelocitiesAsTXT(Parameter* para, int level, int timestep)
    {
        // calc
        int numberNodes = (int)para->getParH(level)->numberOfNodes;
        // write
        UbFileOutputASCII out(para->getFName() + "_VelocitiesASCII_" + std::to_string(level) + "_ID_" +
                              StringUtil::toString<int>(para->getMyProcessID()) + "_t_" + std::to_string(timestep) + ".dat");
        // header
        out.writeString("Level:");
        out.writeInteger(level);
        out.writeLine();
        out.writeString("vX vY vZ");
        out.writeLine();
        // out.writeInteger(numberNodes);
        // out.writeLine();
        for (int index = 0; index < numberNodes; index++) {
            out.writeFloat((float)(para->getParH(level)->velocityX[index] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->velocityY[index] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->velocityZ[index] * para->getVelocityRatio()));
            out.writeLine();
        }
        out.writeLine();
    }
};
#endif
