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
#ifndef CoordNeighborGeoV_H
#define CoordNeighborGeoV_H

#include <vector>
#include <string>
#include <fstream>
#include "Calculation/Calculation.h"


class CoordNeighborGeoV 
{
protected:
    std::ifstream file; 
    unsigned int maxLevel; 
    std::vector<unsigned int> levelSizes;
    std::vector< std::vector<unsigned int> > neighbors;
    std::vector< std::vector< real> > coordinates; 

public:
    CoordNeighborGeoV();
    CoordNeighborGeoV(std::string ad, bool binaer, bool coord);
    ~CoordNeighborGeoV(void);

    void init(bool coord);
    void init_Binary(bool coord);

    unsigned int getLevel();
    unsigned int getSize(unsigned int level);
    std::vector<unsigned int>getVec(unsigned int level);
    void setVec(unsigned int level, std::vector<unsigned int> vec);

    void initalNeighbors(unsigned int *int_ptr, unsigned int level ) const;
    void initalCoords(real *int_ptr, unsigned int level ) const;

protected:
    void skipSpace();
    void readLevelSize(unsigned int level);
    void readNeighbors(unsigned int level);
    void readCoordinates(unsigned int level);
    void resizeVectors();
    void readLevel();

};

#endif
