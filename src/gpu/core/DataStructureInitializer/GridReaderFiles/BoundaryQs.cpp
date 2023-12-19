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
#include "BoundaryQs.h"

#include "Parameter/Parameter.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define QCOLUMNS 27

BoundaryQs::BoundaryQs(std::string path, bool isBinary)
{
    if (isBinary) 
    {
        file.open(path.c_str(), std::ios::in | std::ios::binary);
        checkFileStatus(path);
        init_Binary();
    } 
    else
    {
        file.open(path.c_str(), std::ios::in );
        checkFileStatus(path);
        init();
    }
}

BoundaryQs::BoundaryQs(std::string path, std::shared_ptr<Parameter> para, std::string str, bool isBinary)
{
    if (isBinary) file.open(path.c_str(), std::ios::in | std::ios::binary);
    else file.open(path.c_str(), std::ios::in);

    if (!file) 
        para->setObj(str, false);
    else 
    {
        para->setObj(str, true);
        if(isBinary)
            init_Binary();
        else
            init();
    }
    
}

BoundaryQs::BoundaryQs()
{

}


void BoundaryQs::checkFileStatus(std::string path)
{
    if (!file)
    {
        std::cerr << "can not open q-file: " << path << std::endl;
        exit(1);
    }
}

BoundaryQs::~BoundaryQs()
{
    file.close();
}

void BoundaryQs::init()
{
    std::vector<uint32_t> vec1D_code;
    std::string bufferString;

    file >> maxLevel;
    resizeVectors();

    for (unsigned int level = 0; level <= maxLevel; level++)
    {
        file >> levelSizes[level];
        resizeVectorsPerLevel(level, vec1D_code);
        if (levelSizes[level] == 0)
            continue;
        
        for (unsigned int elem = 0; elem < levelSizes[level]; elem++)
        {
            int columnCounter = 26;
            file >> indices[level][elem];
            file >> vec1D_code[elem];

            while (vec1D_code[elem] != 0)
            {
                if (vec1D_code[elem] % 2 == 1)
                    file >> values[level][columnCounter][elem];
                vec1D_code[elem] /= 2;
                columnCounter--;
            }
            getline(file, bufferString);
        }
        vec1D_code.clear();
    }
}

void BoundaryQs::init_Binary() 
{
    std::vector<uint32_t> vec1D_code;

    std::string bufferString;
    unsigned int bufferInt;
    real bufferDouble;
    uint32_t bufferUint32_t;

    file >> maxLevel;
    resizeVectors();

    for (unsigned int level = 0; level <= maxLevel; level++)
    {
        file >> levelSizes[level];
        resizeVectorsPerLevel(level, vec1D_code);
        if (levelSizes[level] == 0)
            continue;

        for (unsigned int elem = 0; elem < levelSizes[level]; elem++)
        {
            int zaehler = 26;
            file.read((char*)&bufferInt, sizeof(int));
            indices[level][elem] = bufferInt;

            file.read((char*)&bufferUint32_t, sizeof(uint32_t));
            vec1D_code[elem] = bufferUint32_t;
            while (vec1D_code[elem] != 0)
            {
                if (vec1D_code[elem] % 2 == 1)
                {
                    file.read((char*)&bufferDouble, sizeof(double));
                    values[level][zaehler][elem] = bufferDouble;
                }
                vec1D_code[elem] /= 2;
                zaehler--;
            }
            getline(file, bufferString);
        }
        vec1D_code.clear();
    }
}

void BoundaryQs::resizeVectors()
{
    levelSizes.resize(maxLevel + 1);
    values.resize(maxLevel + 1);
    indices.resize(maxLevel + 1);
}

void BoundaryQs::resizeVectorsPerLevel(unsigned int level, std::vector<uint32_t> &vec1D_code)
{
    values[level].resize(QCOLUMNS);
    for (int i = 0; i < QCOLUMNS; i++)
        values[level][i].resize(levelSizes[level], -1);
    indices[level].resize(levelSizes[level]);
    vec1D_code.resize(levelSizes[level]);
}

unsigned int BoundaryQs::getSize(unsigned int level)
{
    return this->levelSizes[level];
}

unsigned int BoundaryQs::getLevel()
{
    return maxLevel;
}


void BoundaryQs::setValuesInVector(std::vector<std::vector<std::vector<real>>> &q27, unsigned int level) const 
{
    for (std::size_t column = 0; column < values[level].size(); column++)
        for (std::size_t index = 0; index < values[level][column].size(); index++)
            q27[level][column].push_back(values[level][column][index]);
}

void BoundaryQs::setValues(real **q27, unsigned int level) const
{
    for (std::size_t column = 0; column < values[level].size(); column++)
        for (std::size_t index = 0; index < values[level][column].size(); index++)
            q27[column][index] = values[level][column][index];
}

void BoundaryQs::setIndexInVector(std::vector<std::vector<int>> &data, unsigned int level) const 
{
    for (std::size_t index = 0; index < indices[level].size(); index++)
        data[level].push_back(indices[level][index]);
}

void BoundaryQs::setIndex(int *data, unsigned int level) const
{
    for (std::size_t index = 0; index < indices[level].size(); index++)
        data[index] = indices[level][index];
}


void BoundaryQs::getQs(std::vector<std::vector<std::vector<real> > > &qs) {
    int j = 0;
    int i = 0;
    for (std::vector<std::vector<std::vector<real> > >::iterator it = values.begin(); it != values.end(); it++) {
        i = 0;
        for (std::vector<std::vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {

            for (std::vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
                qs[j][i].push_back(*it3);
            }
            i++;
        }
        j++;
    }
}

void BoundaryQs::getIndices(std::vector<std::vector<uint> > &indices) {
    for (std::size_t level = 0; level < this->indices.size(); level++)
        for (std::size_t index = 0; index < this->indices[level].size(); index++)
            indices[level].push_back(this->indices[level][index]);
}


//! \}
