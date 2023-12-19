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
#ifndef MATHEMATICA_ASSISTANT_IMP_H
#define MATHEMATICA_ASSISTANT_IMP_H

#include "MathematicaAssistant.h"

class MathematicaFunctionFactory;

class MathematicaAssistantImp : public MathematicaAssistant
{
public:
    virtual void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile) = 0;

protected:
    MathematicaAssistantImp();
    MathematicaAssistantImp(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::string timeStep);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::vector<std::string> timeSteps);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::vector<std::string> timeSteps, std::vector<std::string> basicKernels);
    std::string finalizeListName(std::string basicListName, std::string testName, std::string dataToCalc);

    void addListLogLogPlotToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::string> listNames, std::vector<std::vector<double> > xAxesData, std::vector<std::vector<double> > yAxesData, std::string labelXAxes, std::string labelYAxes);
    void addListOfListsToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::string listName, std::vector<std::vector<double> > listOfLists);

    void addSecondOrderOfAccuracyRef(std::vector<std::vector<double> > &xAxesData, std::vector<std::vector<double> > &yAxesData, std::vector<std::string> &listNames);
    void addFourthOrderOfAccuracyRef(std::vector<std::vector<double> > &xAxesData, std::vector<std::vector<double> > &yAxesData, std::vector<std::string> &listNames);

    std::shared_ptr<MathematicaFunctionFactory> functionFactory;
};
#endif 
//! \}
