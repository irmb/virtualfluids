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
//! \author Martin Schoenherr
//=======================================================================================
#include "ReaderMeasurePoints.h"

#include "Parameter/Parameter.h"
#include "Cuda/CudaMemoryManager.h"

#include <basics/utilities/UbFileInputASCII.h>

using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////
void ReaderMeasurePoints::readMeasurePoints( Parameter* para ) 
{
    UbFileInputASCII in(para->getmeasurePoints());
    int numberOfAllNodes = in.readInteger();
    in.readLine();
    int tempLevel;
    MeasurePoints tempMP;
    //printf("done, init the values...\n");
    for (int u = 0; u < numberOfAllNodes; u++)
    {
        tempMP.name = in.readString();         
        //printf("done, read the name...\n");
        tempMP.k = in.readInteger();
        //printf("done, read k...\n");
        tempLevel = in.readInteger();
        //printf("done, read level...\n");
        in.readLine();
        //printf("done, read the values...\n");
        para->getParH(tempLevel)->MP.push_back(tempMP);
        //printf("done, put it into a vector...\n");
    }
}

////////////////////////////////////////////////////////////////////////////////
void ReaderMeasurePoints::readMeasurePoints(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    // read measure points from file
    ReaderMeasurePoints::readMeasurePoints(para);
    // printf("done, reading the file...\n");
    // level loop
    for (int lev = 0; lev <= para->getMaxLevel(); lev++) {
        // set Memory Size and malloc of the indices and macroscopic values per level
        para->getParH(lev)->numberOfValuesMP = (unsigned int)para->getParH(lev)->MP.size() *
                                               (unsigned int)para->getclockCycleForMP() /
                                               ((unsigned int)para->getTimestepForMP());
        para->getParD(lev)->numberOfValuesMP = para->getParH(lev)->numberOfValuesMP;

        para->getParH(lev)->numberOfPointskMP = (int)para->getParH(lev)->MP.size();
        para->getParD(lev)->numberOfPointskMP = para->getParH(lev)->numberOfPointskMP;

        para->getParH(lev)->memSizeIntkMP = sizeof(unsigned int) * (int)para->getParH(lev)->MP.size();
        para->getParD(lev)->memSizeIntkMP = para->getParH(lev)->memSizeIntkMP;

        para->getParH(lev)->memSizerealkMP = sizeof(real) * para->getParH(lev)->numberOfValuesMP;
        para->getParD(lev)->memSizerealkMP = para->getParH(lev)->memSizerealkMP;

        printf("Level: %d, numberOfValuesMP: %d, memSizeIntkMP: %d, memSizerealkMP: %d\n", lev,
               para->getParH(lev)->numberOfValuesMP, para->getParH(lev)->memSizeIntkMP, para->getParD(lev)->memSizerealkMP);

        cudaMemoryManager->cudaAllocMeasurePointsIndex(lev);

        // loop over all measure points per level
        for (int index = 0; index < (int)para->getParH(lev)->MP.size(); index++) {
            // set indices
            para->getParH(lev)->kMP[index] = para->getParH(lev)->MP[index].k;
        }
        // loop over all measure points per level times MPClockCycle
        for (int index = 0; index < (int)para->getParH(lev)->numberOfValuesMP; index++) {
            // init values
            para->getParH(lev)->VxMP[index] = (real)0.0;
            para->getParH(lev)->VyMP[index] = (real)0.0;
            para->getParH(lev)->VzMP[index] = (real)0.0;
            para->getParH(lev)->RhoMP[index] = (real)0.0;
        }

        // copy indices-arrays
        cudaMemoryManager->cudaCopyMeasurePointsIndex(lev);
    }
}
////////////////////////////////////////////////////////////////////////////////
