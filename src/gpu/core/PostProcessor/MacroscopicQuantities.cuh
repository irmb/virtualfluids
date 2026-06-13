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
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================

#ifndef MacroscopicQuantities_CUH
#define MacroscopicQuantities_CUH

#include <basics/DataTypes.h>

namespace vf::gpu {

void calculateMacroscopicQuantities(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep);

void calculateMacroscopicQuantitiesCompressible(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, const uint* geoD, const uint* neighborX,
                     const uint* neighborY, const uint* neighborZ, unsigned long long numberOfLBnodes,
                     uint numberOfThreads, real* DD, bool isEvenTimestep);

void calculateSubGridScaleFluxesCompressible(const uint* indices, uint numberOfIndices, real* vxvx, real* vxvy, real* vxvz, real* vyvy, real* vyvz, real* vzvz,
                                             real* phix, real* phiy, real* phiz, const uint* geoD, const real* vx,
                                             const real* vy, const real* vz, const real* scalars,
                                             const real* turbulenceViscosities, const real* turbulenceDiffusivities,
                                             const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                             unsigned long long numberOfLBnodes, real* distributions,
                                             real* distributionsScalar, real omega, real omegaDiffusive, uint numberOfThreads,
                                             bool isEvenTimestep, bool computeSubgridScaleFluxesScalar);

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

}

#endif

//! \}
