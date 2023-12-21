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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "MathematicaAssistantFactoryImp.h"

#include "Tests/PhiTest/MathematicaAssistant/PhiMathematicaAssistant.h"
#include "Tests/NyTest/MathematicaAssistant/NyMathematicaAssistant.h"
#include "Tests/L2Norm/MathematicaAssistant/L2NormMathematicaAssistant.h"
#include "Tests/L2NormBetweenKernels/MathematicaAssistant/L2NormBetweenKernelsMathematicaAssistant.h"

#include "Utilities/MathematicaAssistant/TimeAssistant/TimeMathematicaAssistant.h"

std::shared_ptr<MathematicaAssistantFactory> MathematicaAssistantFactoryImp::getNewInstance()
{
    return std::shared_ptr<MathematicaAssistantFactory>(new MathematicaAssistantFactoryImp());
}

std::vector<std::shared_ptr<MathematicaAssistant>> MathematicaAssistantFactoryImp::makeMathematicaAssistants(std::vector<Assistant> types, std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
    std::vector<std::shared_ptr<MathematicaAssistant>> myAssistants;

    for(int i = 0; i < types.size(); i++){
        switch (types.at(i))
        {
        case Phi:
            myAssistants.push_back(PhiMathematicaAssistant::getNewInstance(functionFactory));
            break;
        case Ny:
            myAssistants.push_back(NyMathematicaAssistant::getNewInstance(functionFactory));
            break;
        case L2Norm:
            myAssistants.push_back(L2NormMathematicaAssistant::getNewInstance(functionFactory));
            break;
        case L2NormBetweenKernels:
            myAssistants.push_back(L2NormBetweenKernelsMathematicaAssistant::getNewInstance(functionFactory));
            break;
        case Time:
            myAssistants.push_back(TimeMathematicaAssistant::getNewInstance(functionFactory));
            break;
        default:
            break;
        }

    }


    return myAssistants;
}

MathematicaAssistantFactoryImp::MathematicaAssistantFactoryImp()
{
}




//! \}
