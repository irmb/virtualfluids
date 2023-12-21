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
#include "MathematicaAssistantImp.h"

#include "Utilities/DataPoint/DataPoint.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactory.h"

#include <cmath>


MathematicaAssistantImp::MathematicaAssistantImp()
{
}

MathematicaAssistantImp::MathematicaAssistantImp(std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
    this->functionFactory = functionFactory;
}

std::vector<std::string> MathematicaAssistantImp::finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc)
{
    std::vector<std::string> finalListNames;

    for (int i = 0; i < basicListNames.size(); i++)
        finalListNames.push_back(finalizeListName(basicListNames.at(i), testName, dataToCalc));

    return finalListNames;
}

std::vector<std::string> MathematicaAssistantImp::finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData)
{
    std::vector<std::string> finalListNames;

    for (int i = 0; i < basicListNames.size(); i++)
        finalListNames.push_back(testName + basicListNames.at(i) + dataToCalc + normalizeData);

    return finalListNames;
}

std::vector<std::string> MathematicaAssistantImp::finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::string timeStep)
{
    std::vector<std::string> finalListNames;

    for (int i = 0; i < basicListNames.size(); i++)
        finalListNames.push_back(testName + basicListNames.at(i) + dataToCalc + normalizeData + timeStep);

    return finalListNames;
}

std::vector<std::string> MathematicaAssistantImp::finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::vector<std::string> timeSteps)
{
    std::vector<std::string> finalListNames;

    for (int i = 0; i < basicListNames.size(); i++)
        finalListNames.push_back(testName + basicListNames.at(i) + dataToCalc + normalizeData + timeSteps.at(i));

    return finalListNames;
}

std::vector<std::string> MathematicaAssistantImp::finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::vector<std::string> timeSteps, std::vector<std::string> basicKernels)
{
    std::vector<std::string> finalListNames;

    for (int i = 0; i < basicListNames.size(); i++)
        finalListNames.push_back(testName + basicListNames.at(i) + dataToCalc + normalizeData + timeSteps.at(i) + basicKernels.at(i));

    return finalListNames;
}

std::string MathematicaAssistantImp::finalizeListName(std::string basicListName, std::string testName, std::string dataToCalc)
{
    return testName + basicListName + dataToCalc;
}

void MathematicaAssistantImp::addListLogLogPlotToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::string> listNames, std::vector<std::vector<double>> xAxesData, std::vector<std::vector<double>> yAxesData, std::string labelXAxes, std::string labelYAxes)
{
    std::vector<std::vector<std::shared_ptr<DataPoint> > > dataPointGroup;

    for (int i = 0; i < xAxesData.size(); i++) {
        std::vector<std::shared_ptr<DataPoint> > dataPoints;
        for (int j = 0; j < xAxesData.at(i).size(); j++)
            dataPoints.push_back(DataPoint::getNewInstance(xAxesData.at(i).at(j), yAxesData.at(i).at(j)));
        dataPointGroup.push_back(dataPoints);
    }
    std::vector<std::shared_ptr<MathematicaPointList> > pointList;
    for (int i = 0; i < dataPointGroup.size(); i++) {
        std::shared_ptr<MathematicaPointList> aPointList = functionFactory->makeMathematicaPointList(aMathmaticaFile, listNames.at(i), dataPointGroup.at(i));
        pointList.push_back(aPointList);
    }
    std::shared_ptr<MathematicaListPlot> listLogLogPlot = functionFactory->makeMathematicaListPlot(aMathmaticaFile, pointList, "ListLogLogPlot", labelXAxes, labelYAxes);
}

void MathematicaAssistantImp::addListOfListsToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::string listNames, std::vector<std::vector<double>> listOfLists)
{
    functionFactory->makeMathematicaListOfLists(aMathmaticaFile, listNames, listOfLists);
}

void MathematicaAssistantImp::addSecondOrderOfAccuracyRef(std::vector<std::vector<double>>& xAxesData, std::vector<std::vector<double>>& yAxesData, std::vector<std::string>& listNames)
{
    int  maxLength = 0;
    int maxLengthAtNumber = 0;
    for (int i = 0; i < xAxesData.size(); i++) {
        if (xAxesData.at(i).size() > maxLength) {
            maxLength = xAxesData.at(i).size();
            maxLengthAtNumber = i;
        }
    }
    std::vector<double> x = xAxesData.at(maxLengthAtNumber);
    
    double maxData = 0.0;
    for (int i = 0; i < yAxesData.size(); i++) {
        for (int j = 0; j < yAxesData.at(i).size(); j++) {
            if (yAxesData.at(i).at(j) > maxData)
                maxData = yAxesData.at(i).at(j);
        }
    }

    std::vector<double> sec = { maxData * 10.0 };
    for (int l = 1; l < x.size(); l++) 
        sec.push_back(sec.at(l - 1) / exp(-2.0 * log(x.at(l - 1) / x.at(l))));
    xAxesData.push_back(x);
    yAxesData.push_back(sec);
    listNames.push_back("SecondOrderOfAccuracy");
}

void MathematicaAssistantImp::addFourthOrderOfAccuracyRef(std::vector<std::vector<double>>& xAxesData, std::vector<std::vector<double>>& yAxesData, std::vector<std::string>& listNames)
{
    int  maxLength = 0;
    int maxLengthAtNumber = 0;
    for (int i = 0; i < xAxesData.size(); i++) {
        if (xAxesData.at(i).size() > maxLength) {
            maxLength = xAxesData.at(i).size();
            maxLengthAtNumber = i;
        }
    }
    std::vector<double> x = xAxesData.at(maxLengthAtNumber);

    double minData = yAxesData.at(0).at(0);
    for (int i = 1; i < yAxesData.size(); i++) {
            if (yAxesData.at(i).at(0) < minData)
                minData = yAxesData.at(i).at(0);
    }
    std::vector<double> fourth = { minData / 10.0 };
    for (int l = 1; l < x.size(); l++) 
        fourth.push_back(fourth.at(l - 1) / std::exp(-4.0 * std::log(x.at(l - 1) / x.at(l))));


    xAxesData.push_back(x);
    yAxesData.push_back(fourth);
    listNames.push_back("FourthOrderOfAccuracy");
}

//! \}
