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
#include "MathematicaPointList.h"

#include "Utilities/DataPoint/DataPoint.h"

std::shared_ptr<MathematicaPointList> MathematicaPointList::getNewInstance(std::string listName, std::vector< std::shared_ptr< DataPoint>> plotData)
{
    return std::shared_ptr<MathematicaPointList>(new MathematicaPointList(listName, plotData));
}

std::string MathematicaPointList::getListName()
{
    return listName;
}

MathematicaPointList::MathematicaPointList(std::string listName, std::vector< std::shared_ptr< DataPoint>> plotData) : listName(listName)
{
    mathematicaFunction << listName << "= {";
    for (int i = 0; i < plotData.size(); i++) {
        if (i > 0)
            mathematicaFunction << ", ";
        mathematicaFunction << "{" << plotData.at(i)->getX() << ", " << plotData.at(i)->getY() << "}";
    }

    mathematicaFunction << "}";
}

MathematicaPointList::MathematicaPointList()
{
}



//! \}
