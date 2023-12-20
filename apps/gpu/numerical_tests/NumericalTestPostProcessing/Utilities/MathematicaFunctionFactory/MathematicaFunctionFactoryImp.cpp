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
#include "MathematicaFunctionFactoryImp.h"

#include "Utilities/MathematicaFile/MathematicaFile.h"

#include "Utilities/MathematicaFunction/LinePlot/MathematicaListPlotImp.h"
#include "Utilities/MathematicaFunction/PointList/MathematicaPointListImp.h"
#include "Utilities/MathematicaFunction/ListOfLists/MathematicaListOfListsImp.h"

std::shared_ptr<MathematicaFunctionFactory> MathematicaFunctionFactoryImp::getNewInstance()
{
    return std::shared_ptr<MathematicaFunctionFactory>(new MathematicaFunctionFactoryImp());
}

std::shared_ptr<MathematicaPointList> MathematicaFunctionFactoryImp::makeMathematicaPointList(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::shared_ptr<DataPoint>> plotData)
{
    std::shared_ptr<MathematicaPointList> mathPointList = MathematicaPointListImp::getNewInstance(listName, plotData);
    file->addMathematicaFunction(mathPointList);
    return mathPointList;
}

std::shared_ptr<MathematicaListPlot> MathematicaFunctionFactoryImp::makeMathematicaListPlot(std::shared_ptr<MathematicaFile> file, std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes)
{
    std::shared_ptr<MathematicaListPlot> listLinePlot = MathematicaListPlotImp::getNewInstance(pointList, plotType, xAxes, yAxes);
    file->addMathematicaFunction(listLinePlot);
    return listLinePlot;
}

std::shared_ptr<MathematicaListOfLists> MathematicaFunctionFactoryImp::makeMathematicaListOfLists(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::vector<double>> listOfLists)
{
    std::shared_ptr<MathematicaListOfLists> listOfList = MathematicaListOfListsImp::getNewInstance(listName, listOfLists);
    file->addMathematicaFunction(listOfList);
    return listOfList;
}

MathematicaFunctionFactoryImp::MathematicaFunctionFactoryImp()
{
}

//! \}
