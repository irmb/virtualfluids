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
#include "MeasuredPoints.h"
#include <cstdlib>
#include <iostream>

//using namespace std;

MeasuredPoints::MeasuredPoints(void)
{
}

MeasuredPoints::MeasuredPoints(std::string ad){
    file.open(ad.c_str(), std::ios::in | std::ios::binary);

    if (!file) {
        std::cerr << "Fehler beim Oeffnen Measured Points" << std::endl;
            exit(1);
    }

    this->init();        

}

MeasuredPoints::~MeasuredPoints(void)
{
}




void MeasuredPoints::init() {
    
    std::string bufferString;
    unsigned int bufferInt;

    getline(file,bufferString);
    readLevel();

    this->levelSizes.resize(maxLevel);
    this->points.resize(maxLevel);

    for (uint i=0; i<maxLevel; i++) {
        getline(file,bufferString);
        bufferInt = atoi(bufferString.c_str()); 

        this->levelSizes[i]=bufferInt;

        this->points[i].resize(levelSizes[i]);
        if(levelSizes[i] != 0) {
            for ( uint j = 0; j < levelSizes[i]; j++) {
                getline(file,bufferString);
                bufferInt = atoi(bufferString.c_str()); 
                this->points[i][j]=bufferInt;
            }
        }


    }



}