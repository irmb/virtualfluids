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
#include "MathematicaListPlotImp.h"

#include "Utilities/MathematicaFunction/PointList/MathematicaPointList.h"

std::shared_ptr<MathematicaListPlot> MathematicaListPlotImp::getNewInstance(std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes)
{
    return std::shared_ptr<MathematicaListPlotImp>(new MathematicaListPlotImp(pointList, plotType, xAxes, yAxes));
}

MathematicaListPlotImp::MathematicaListPlotImp(std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes)
{
    mathematicaFunction << plotType << "[{";
    for (int i = 0; i < pointList.size(); i++) {
        if(i < pointList.size() - 1)
            mathematicaFunction << pointList.at(i)->getListName() << ", ";
        else
            mathematicaFunction << pointList.at(i)->getListName() << "}";
    }
    mathematicaFunction << ", PlotLegends -> {\"";
    for (int i = 0; i < pointList.size(); i++) {
        if (i < pointList.size() - 1)
            mathematicaFunction << pointList.at(i)->getListName() << "\", \"";
        else
            mathematicaFunction << pointList.at(i)->getListName() << "\"}";
    }
    mathematicaFunction << ", AxesLabel -> {\"" << xAxes << "\", \"" << yAxes << "\"}, Joined -> True, PlotMarkers->Automatic, PlotStyle -> Dashed]";
}

MathematicaListPlotImp::MathematicaListPlotImp()
{
}
//! \}
