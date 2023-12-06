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
        para->getParH(tempLevel)->MeasurePointVector.push_back(tempMP);
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
        para->getParH(lev)->numberOfMeasurePoints = (unsigned int)para->getParH(lev)->MeasurePointVector.size() *
                                               (unsigned int)para->getclockCycleForMeasurePoints() /
                                               ((unsigned int)para->getTimestepForMeasurePoints());
        para->getParD(lev)->numberOfMeasurePoints = para->getParH(lev)->numberOfMeasurePoints;

        para->getParH(lev)->numberOfMeasurePoints = (int)para->getParH(lev)->MeasurePointVector.size();
        para->getParD(lev)->numberOfMeasurePoints = para->getParH(lev)->numberOfMeasurePoints;

        para->getParH(lev)->memSizeIntegerMeasurePoints = sizeof(unsigned int) * (int)para->getParH(lev)->MeasurePointVector.size();
        para->getParD(lev)->memSizeIntegerMeasurePoints = para->getParH(lev)->memSizeIntegerMeasurePoints;

        para->getParH(lev)->memSizeRealMeasurePoints = sizeof(real) * para->getParH(lev)->numberOfMeasurePoints;
        para->getParD(lev)->memSizeRealMeasurePoints = para->getParH(lev)->memSizeRealMeasurePoints;

        printf("Level: %d, numberOfValuesMP: %d, memSizeIntkMP: %d, memSizerealkMP: %d\n", lev,
               para->getParH(lev)->numberOfMeasurePoints, para->getParH(lev)->memSizeIntegerMeasurePoints, para->getParD(lev)->memSizeRealMeasurePoints);

        cudaMemoryManager->cudaAllocMeasurePointsIndex(lev);

        // loop over all measure points per level
        for (int index = 0; index < (int)para->getParH(lev)->MeasurePointVector.size(); index++) {
            // set indices
            para->getParH(lev)->indicesOfMeasurePoints[index] = para->getParH(lev)->MeasurePointVector[index].k;
        }
        // loop over all measure points per level times MPClockCycle
        for (int index = 0; index < (int)para->getParH(lev)->numberOfMeasurePoints; index++) {
            // init values
            para->getParH(lev)->velocityInXdirectionAtMeasurePoints[index] = (real)0.0;
            para->getParH(lev)->velocityInYdirectionAtMeasurePoints[index] = (real)0.0;
            para->getParH(lev)->velocityInZdirectionAtMeasurePoints[index] = (real)0.0;
            para->getParH(lev)->densityAtMeasurePoints[index] = (real)0.0;
        }

        // copy indices-arrays
        cudaMemoryManager->cudaCopyMeasurePointsIndex(lev);
    }
}
////////////////////////////////////////////////////////////////////////////////
