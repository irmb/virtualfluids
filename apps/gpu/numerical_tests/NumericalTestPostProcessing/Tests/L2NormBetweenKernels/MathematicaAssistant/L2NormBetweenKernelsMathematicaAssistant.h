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
#ifndef L2NORM_BK_MATHEMATICA_ASSISTANT_H
#define L2NORM_BK_MATHEMATICA_ASSISTANT_H

#include "Utilities/MathematicaAssistant/MathematicaAssistantImp.h"

class MathematicaFunctionFactory;
class L2NormBetweenKernelsLogFileData;

struct SortedDataL2NormBetweenKernels {
    std::vector<std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > > testLogFileData;
    std::vector<std::vector<std::string> > basicListNames;
};

class L2NormBetweenKernelsMathematicaAssistant : public MathematicaAssistantImp
{
public:
    static std::shared_ptr<L2NormBetweenKernelsMathematicaAssistant> getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

    
private:
    L2NormBetweenKernelsMathematicaAssistant();
    L2NormBetweenKernelsMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    bool checkTestParameter(std::shared_ptr<L2NormBetweenKernelsLogFileData> logFileData1, std::shared_ptr<L2NormBetweenKernelsLogFileData> logFileData2);
    std::shared_ptr<SortedDataL2NormBetweenKernels> sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData);

    void makeL2NormMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2NormBetweenKernels> sortedData);
    void makeL2NormBetweenKernelsMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2NormBetweenKernels> sortedData);

    
};
#endif 
//! \}
