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
//! \addtogroup gpu_geometries geometries
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "Conglomerate.h"

#include <memory>

SPtr<Conglomerate> Conglomerate::makeShared()
{
    return std::make_shared<Conglomerate>();
}

void Conglomerate::add(SPtr<Object> object)
{
    if (numberOfAddObjects < MAX_NUMBER_OF_OBJECTS)
    {
        addObjects[numberOfAddObjects] = object;
        numberOfAddObjects++;
    } else
        printf("[WARNING] max numbers of %d reached! Object was not added.\n", MAX_NUMBER_OF_OBJECTS);
}

void Conglomerate::subtract(SPtr<Object> object)
{
    if (numberOfSubtractObjects < MAX_NUMBER_OF_OBJECTS)
    {
        subtractObjects[numberOfSubtractObjects] = object;
        numberOfSubtractObjects++;
    }
    else
        printf("[WARNING] max numbers of %d reached! Object was not added.\n", MAX_NUMBER_OF_OBJECTS);
}

SPtr<Object> Conglomerate::clone() const
{
    auto conglomerate = std::make_shared<Conglomerate>();
    for (uint i = 0; i < numberOfAddObjects; i++)
        conglomerate->add(addObjects[i]->clone());

    for (uint i = 0; i < numberOfSubtractObjects; i++)
        conglomerate->subtract(subtractObjects[i]->clone());

    return conglomerate;
}

double Conglomerate::getX1Centroid() const
{
    return (getX1Minimum() + getX1Maximum()) * 0.5;
}

double Conglomerate::getX1Minimum() const
{
    double minimum = addObjects[0]->getX1Minimum();
    for(uint i = 1; i < numberOfAddObjects; i++)
        minimum = getMinimum(minimum, addObjects[i]->getX1Minimum());
    return minimum;
}

double Conglomerate::getX1Maximum() const
{
    double maximum = addObjects[0]->getX1Maximum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        maximum = getMaximum(maximum, addObjects[i]->getX1Maximum());
    return maximum;
}

double Conglomerate::getX2Centroid() const
{
    return (getX2Minimum() + getX2Maximum()) * 0.5;
}

double Conglomerate::getX2Minimum() const
{
    double minimum = addObjects[0]->getX2Minimum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        minimum = getMinimum(minimum, addObjects[i]->getX2Minimum());
    return minimum;
}    

double Conglomerate::getX2Maximum() const
{
    double maximum = addObjects[0]->getX2Maximum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        maximum = getMaximum(maximum, addObjects[i]->getX2Maximum());
    return maximum;
}

double Conglomerate::getX3Centroid() const
{
    return (getX3Minimum() + getX3Maximum()) * 0.5;
}

double Conglomerate::getX3Minimum() const
{
    double minimum = addObjects[0]->getX3Minimum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        minimum = getMinimum(minimum, addObjects[i]->getX3Minimum());
    return minimum;
}    

double Conglomerate::getX3Maximum() const
{
    double maximum = addObjects[0]->getX3Maximum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        maximum = getMaximum(maximum, addObjects[i]->getX3Maximum());
    return maximum;
}

double Conglomerate::getMinimum(double val1, double val2)
{
    if (val1 > val2)
        return val2;
    return val1;
}

double Conglomerate::getMaximum(double val1, double val2)
{
    if (val1 < val2)
        return val2;
    return val1;
}

bool Conglomerate::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
    for (uint i = 0; i < numberOfSubtractObjects; i++)
        if (subtractObjects[i]->isPointInObject(x1, x2, x3, minOffset, maxOffset))
            return false;

    for (uint i = 0; i < numberOfAddObjects; i++)
        if (addObjects[i]->isPointInObject(x1, x2, x3, minOffset, maxOffset))
            return true;

    return false;
}

void Conglomerate::changeSizeByDelta(double delta)
{
    for (uint i = 0; i < numberOfAddObjects; i++)
        addObjects[i]->changeSizeByDelta(delta);

    for (uint i = 0; i < numberOfSubtractObjects; i++)
        subtractObjects[i]->changeSizeByDelta(delta);
}

void Conglomerate::findInnerNodes(SPtr<GridImp> grid)
{
    for (uint i = 0; i < numberOfAddObjects; i++)
        addObjects[i]->findInnerNodes(grid);

    if( numberOfSubtractObjects > 0 )
        VF_LOG_WARNING("Warning: Conglomerate::substract() is currently nut fully implemented!");
}

//! \}
