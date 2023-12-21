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
#include "MathematicaListOfListsImp.h"

#include <iomanip>
#include <limits>

std::shared_ptr<MathematicaListOfLists> MathematicaListOfListsImp::getNewInstance(std::string listName, std::vector<std::vector<double>> listOfLists)
{
    return std::shared_ptr<MathematicaListOfLists>(new MathematicaListOfListsImp(listName, listOfLists));
}

std::string MathematicaListOfListsImp::getListName()
{
    return listName;
}

MathematicaListOfListsImp::MathematicaListOfListsImp(std::string listName, std::vector<std::vector<double>> listOfLists) : listName(listName), listOfLists(listOfLists)
{
    mathematicaFunction << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    mathematicaFunction << listName << "= {";

    for (int i = 0; i < listOfLists.size(); i++) {
        mathematicaFunction << "{";
        for (int j = 0; j < listOfLists.at(i).size(); j++){
            if (j > 0)
                mathematicaFunction << ", ";
            mathematicaFunction << listOfLists.at(i).at(j);
        }
        mathematicaFunction << "}";
        if(i < listOfLists.size() - 1)
            mathematicaFunction << ", ";
    }
    mathematicaFunction << "};";
}

MathematicaListOfListsImp::MathematicaListOfListsImp()
{
}

//! \}
