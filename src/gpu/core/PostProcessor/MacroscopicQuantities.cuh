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
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================

#ifndef MacroscopicQuantities_CUH
#define MacroscopicQuantities_CUH

#include <basics/DataTypes.h>

void calculateMacroscopicQuantities(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep);

void calculateMacroscopicQuantitiesCompressible(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                     unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                     unsigned int numberOfThreads, real* DD, bool isEvenTimestep);

void calculateMean(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep);

void calculateMeanCompressible(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                     unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                     unsigned int numberOfThreads, real* DD, bool isEvenTimestep);

void calculateMeanCompressibleAdvectionDiffusion(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, real* concD, unsigned int* geoD,
                     unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                     unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD, real* DD_AD,
                     bool isEvenTimestep);

void calculateMacrosopicMean(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                    unsigned int* neighborY, unsigned int* neighborZ, unsigned int tdiff, unsigned long long numberOfLBnodes,
                    unsigned int numberOfThreads, bool isEvenTimestep);

void ResetMeanValuesSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned long long numberOfLBnodes,
                         unsigned int numberOfThreads, bool isEvenTimestep);

void ResetMeanValuesAD27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, real* concD,
                         unsigned long long numberOfLBnodes, unsigned int numberOfThreads, bool isEvenTimestep);

void calculateMeasurePoints(real* vxMP, real* vyMP, real* vzMP, real* rhoMP, unsigned int* kMP,
                           unsigned int numberOfPointskMP, unsigned int MPClockCycle, unsigned int t, unsigned int* geoD,
                           unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                           unsigned long long numberOfLBnodes, real* DD, unsigned int numberOfThreads, bool isEvenTimestep);

#endif
