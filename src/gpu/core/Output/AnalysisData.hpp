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
#ifndef ANALYSIS_DATA_H
#define ANALYSIS_DATA_H

#include <basics/StringUtilities/StringUtil.h>
#include <basics/constants/NumericConstants.h>
#include <basics/utilities/UbFileOutputASCII.h>

#include "Parameter/Parameter.h"

class AnalysisData
{
public:
    static void writeAnalysisDataXSP(Parameter* para, unsigned int timestep)
    {
        UbFileOutputASCII out(para->getFName() + "_AD_X_" + StringUtil::toString<int>(timestep) + ".dat");

        real level = 0; // uniform
        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = vf::basics::constant::c1o1 / pow(2, level);
        real halfDx = deltaX / vf::basics::constant::c2o1;
        real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / vf::basics::constant::c2o1;

        for (int u = 0; u < numberNodes; u++) {
            if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
                ((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
                ((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u])) {
                out.writeDouble((double)(para->getParH(level)->velocityX[u]));
            }
        }
    }

    static void writeAnalysisDataYSP(Parameter* para, unsigned int timestep)
    {
        UbFileOutputASCII out(para->getFName() + "_AD_Y_" + StringUtil::toString<int>(timestep) + ".dat");

        real level = 0; // uniform
        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = vf::basics::constant::c1o1 / pow(2, level);
        real halfDx = deltaX / vf::basics::constant::c2o1;
        real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / vf::basics::constant::c2o1;

        for (int u = 0; u < numberNodes; u++) {
            if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
                ((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
                ((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u])) {
                out.writeDouble((double)(para->getParH(level)->velocityY[u]));
            }
        }
    }

    static void writeAnalysisDataZSP(Parameter* para, unsigned int timestep)
    {
        UbFileOutputASCII out(para->getFName() + "_AD_Z_" + StringUtil::toString<int>(timestep) + ".dat");

        real level = 0; // uniform
        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = vf::basics::constant::c1o1 / pow(2, level);
        real halfDx = deltaX / vf::basics::constant::c2o1;
        real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / vf::basics::constant::c2o1;

        for (int u = 0; u < numberNodes; u++) {
            if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
                ((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
                ((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u])) {
                out.writeDouble((double)(para->getParH(level)->velocityZ[u]));
            }
        }
    }
};
#endif

//! \}
