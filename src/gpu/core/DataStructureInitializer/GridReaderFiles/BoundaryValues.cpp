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
//! \addtogroup gpu_DataStructureInitializer DataStructureInitializer
//! \ingroup gpu_core core
//! \{
//=======================================================================================
#include "BoundaryValues.h"

#include "Parameter/Parameter.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

BoundaryValues::BoundaryValues(std::string path)
{
    file.open(path.c_str(), std::ios::in);

    if (!file) {
        std::cerr << "error can not open value file:" << path << std::endl;
        exit(1);
    }
    init();
}

BoundaryValues::BoundaryValues(std::string path, std::shared_ptr<Parameter> para, std::string str)
{
    if (!file) {
        std::cerr << "error can not open value file:" << path << std::endl;
        exit(1);
    }

    file.open(path.c_str(), std::ios::in);
    if (!file) {
        para->setObj(str, false);
    } else {
        init();
        para->setObj(str, true);
    }
}


BoundaryValues::BoundaryValues(int neighbor, std::shared_ptr<Parameter> para, std::string sor, std::string direction)
{
    if (direction=="X")
    {
        std::string ad = para->getPossNeighborFilesX(sor)[neighbor];
        file.open(ad.c_str(), std::ios::in);

        if (file.fail()) {
            para->setIsNeighborX(false);
        } else {
            para->setIsNeighborX(true);
            init();
        }
    } 
    else if (direction=="Y")
    {
        std::string ad = para->getPossNeighborFilesY(sor)[neighbor];
        file.open(ad.c_str(), std::ios::in);

        if (file.fail()) {
            para->setIsNeighborY(false);
        } else {
            para->setIsNeighborY(true);
            init();
        }
    }
    else
    {
        std::string ad = para->getPossNeighborFilesZ(sor)[neighbor];
        file.open(ad.c_str(), std::ios::in);

        if (file.fail()) {
            para->setIsNeighborZ(false);
        } else {
            para->setIsNeighborZ(true);
            init();
        }
    }
}


BoundaryValues::~BoundaryValues(void)
{
    file.close();
}


int BoundaryValues::getNumberOfColumns()
{
    if (boundaryCondition == "velocity")
        return 3;
    if (boundaryCondition == "noSlip")
        return 3;
    if (boundaryCondition == "pressure")
        return 2;
    if (boundaryCondition == "processor")
        return 0;
    if (boundaryCondition == "concentration")
        return 0;
    else
        return -1;
}


void BoundaryValues::init() 
{
    readBC();
    readNumberOfLevels();
    resizeVectors();

    int maxColumn = getNumberOfColumns();

    if (maxColumn == -1)
        initalVectorsWithSingleZero();
    else
        initalVectors(maxColumn);
}

void BoundaryValues::initalVectors(unsigned int maxColumn)
{
    for (unsigned int level = 0; level <= maxLevel; level++)
    {
        readLevelSize(level);

        if (levelSizes[level] == 0)
        {
            skipLine();
            continue;
        }

        resizeVectorsPerLevel(level, maxColumn);

        for (unsigned int index = 0; index < levelSizes[level]; index++)
            readData(level, index, maxColumn);

    }
}

void BoundaryValues::readData(unsigned int level, int index, unsigned int maxColumn)
{
    file >> indices[level][index];
    for (unsigned int column = 0; column < maxColumn; column++)
        file >> values[level][column][index];
}

void BoundaryValues::resizeVectorsPerLevel(unsigned int level, unsigned int maxColumn)
{
    values[level].resize(maxColumn);
    indices[level].resize(levelSizes[level]);
    for (unsigned int column = 0; column < maxColumn; column++)
        values[level][column].resize(levelSizes[level]);
}

void BoundaryValues::skipLine()
{
    std::string bufferString;
    getline(file, bufferString);
}

void BoundaryValues::readLevelSize(unsigned int level)
{
    file >> levelSizes[level];
}

void BoundaryValues::initalVectorsWithSingleZero()
{
    indices = { {0} };
    values = { { {0.0} } };
}

void BoundaryValues::readNumberOfLevels()
{
    file >> maxLevel;
}

void BoundaryValues::readBC()
{
    file >> boundaryCondition;
}

void BoundaryValues::resizeVectors()
{
    levelSizes.resize(maxLevel + 1);
    values.resize(maxLevel + 1);
    indices.resize(maxLevel + 1);
}

void BoundaryValues::setBoundarys(std::vector<std::vector<std::vector<real> > > &qs) const
{
    for (unsigned int level = 0; level < values.size(); level++)
        for (unsigned int index = 0; index < values[level].size(); index++)
            for (unsigned int value = 0; value < values[level][index].size(); value++)
                qs[level][index].push_back(values[level][index][value]);
}

void BoundaryValues::setValues(real* velo, unsigned int level, unsigned int column) const
{
    for (std::size_t index = 0; index < values[level][column].size(); index++)
        velo[index] = values[level][column][index];
}

void BoundaryValues::initIndex(/*unsigned*/ int *ptr, unsigned int level)
{
    for (std::size_t i = 0; i < indices[level].size(); i++)
        ptr[i] = indices[level][i];
}

unsigned int BoundaryValues::getLevel()
{
    return maxLevel;
}

unsigned int BoundaryValues::getSize(unsigned int level)
{
    return this->levelSizes[level];
}

std::string BoundaryValues::getBoundaryCondition() 
{
    return this->boundaryCondition;
}


void BoundaryValues::setProcNeighbor(bool pN)
{
    procNeighbor = pN;
}

bool BoundaryValues::getProcNeighbor()
{
    return procNeighbor;
}


void BoundaryValues::setPressValues(real *RhoBC, int* kN, int level) const
{
    for (std::size_t column = 0; column < values[level].size(); column++) {
        for (std::size_t index = 0; index < values[level][column].size(); index++) {
            if (column == 0) RhoBC[index] = values[level][column][index];
            if (column == 1) kN[index] = (int)values[level][column][index];
        }
    }
}

void BoundaryValues::setVelocityValues(real *vx, real *vy, real *vz, int level) const
{
    for (std::size_t column = 0; column < values[level].size(); column++) {
        for (std::size_t index = 0; index < values[level][column].size(); index++) {
            if (column == 0) vx[index] = values[level][column][index];
            if (column == 1) vy[index] = values[level][column][index];
            if (column == 2) vz[index] = values[level][column][index];
        }
    }
}

void BoundaryValues::setOutflowValues(real *RhoBC, int* kN, int level) const
{
    for (std::size_t column = 0; column < values[level].size(); column++) {
        for (std::size_t index = 0; index < values[level][column].size(); index++) {
            if (column == 0) RhoBC[index] = values[level][column][index];
            if (column == 1) kN[index] = (int)values[level][column][index];
        }
    }
}


//! \}
