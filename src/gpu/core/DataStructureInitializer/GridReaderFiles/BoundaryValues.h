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
#ifndef BoundaryValues_H
#define BoundaryValues_H

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "Calculation/Calculation.h"

class Parameter;

class BoundaryValues
{
private:
    std::string boundaryCondition;
    bool procNeighbor;

    std::vector< std::vector<std::vector<real> > >values;
    std::vector< std::vector<unsigned int> >indices;
    std::vector<unsigned int> levelSizes;

    std::ifstream file;
    unsigned int maxLevel;
    
public:
    BoundaryValues(std::string path);
    BoundaryValues(std::string path, std::shared_ptr<Parameter> para, std::string str);
    BoundaryValues(int neighbor, std::shared_ptr<Parameter> para, std::string sor, std::string dir);
    ~BoundaryValues();

    unsigned int getLevel();
    unsigned int getSize(unsigned int level);
    std::string getBoundaryCondition();
    void setBoundarys(std::vector<std::vector<std::vector<real> > > &qs) const;
    void setValues(real* velo, unsigned int level, unsigned int column) const;
    void initIndex(/*unsigned*/ int *ptr, unsigned int level);

    void setProcNeighbor(bool pN);
    bool getProcNeighbor();

    void setPressValues(real *RhoBC, int* kN, int level) const;
    void setVelocityValues(real *vx, real *vy, real *vz, int level) const;
    void setOutflowValues(real *RhoBC, int* kN, int level) const;

private:
    void init();
    int getNumberOfColumns();
    void initalVectors(unsigned int maxColumn);
    void readData(unsigned int level, int index, unsigned int maxColumn);
    void resizeVectorsPerLevel(unsigned int level, unsigned int maxColumn);
    void skipLine();
    void readLevelSize(unsigned int level);
    void initalVectorsWithSingleZero();
    void readNumberOfLevels();
    void readBC();
    void resizeVectors();


};

#endif

//! \}
