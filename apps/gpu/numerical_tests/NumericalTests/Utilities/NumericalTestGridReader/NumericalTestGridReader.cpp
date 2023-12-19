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
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "NumericalTestGridReader.h"

#include "gpu/core/Parameter/Parameter.h"

#include "Utilities/InitialCondition/InitialCondition.h"

#define _USE_MATH_DEFINES
#include <math.h>


std::shared_ptr<NumericalTestGridReader> NumericalTestGridReader::getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager)
{
    return std::shared_ptr<NumericalTestGridReader>(new NumericalTestGridReader(para, initialCondition, cudaManager));
}

void NumericalTestGridReader::setInitialNodeValues(uint numberOfNodes, int level) const
{
    initialCondition->init(level);
    for (uint j = 0; j <= numberOfNodes; j++){
        para->getParH(level)->velocityX[j] = initialCondition->getInitVX(j, level);
        para->getParH(level)->velocityY[j] = initialCondition->getInitVY(j, level);
        para->getParH(level)->velocityZ[j] = initialCondition->getInitVZ(j, level);
        para->getParH(level)->rho[j] = initialCondition->getInitROH(j, level);
        para->getParH(level)->pressure[j] = initialCondition->getInitPRESS(j, level);
    }
}

NumericalTestGridReader::NumericalTestGridReader(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager) : GridReader(FILEFORMAT::BINARY, para, cudaManager), initialCondition(initialCondition)
{

}
//! \}
