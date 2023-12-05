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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//=======================================================================================
#ifndef BoundaryQs_H
#define BoundaryQs_H

#include <vector>
#include <string>
#include <fstream>
#include <memory>

#include "Calculation/Calculation.h"

class Parameter;

class BoundaryQs
{
public:
    BoundaryQs();
    BoundaryQs(std::string q, bool binaer);
    BoundaryQs(std::string q, std::shared_ptr<Parameter> para, std::string str, bool binaer);
    ~BoundaryQs(void);

public:
    unsigned int getSize(unsigned int level);
    unsigned int getLevel();


private:    
    void checkFileStatus(std::string path);
    void init();
    void resizeVectors();
    void resizeVectorsPerLevel(unsigned int level, std::vector<uint32_t> &vec1D_code);

    void init_Binary();

public:
    void setIndexInVector(std::vector<std::vector<int>> &data, unsigned int level) const;
    void setValuesInVector(std::vector<std::vector<std::vector<real>>> &q27, unsigned int level) const;
    void setIndex(int *indices, unsigned int level) const;
    void setValues(real** q27, unsigned int level) const;
    void getQs(std::vector<std::vector<std::vector<real> > > &qs);
    void getIndices(std::vector<std::vector<uint> > &indices);
    //void initArray(real* ptr, unsigned int level, unsigned int column);
    //void initIndex(int *ptr, unsigned int level);

private:
    std::vector< std::vector<std::vector<real> > >values;
    std::vector< std::vector<unsigned int> >indices;

    std::ifstream file;
    unsigned int maxLevel;
    std::vector<unsigned int> levelSizes;

};

#endif