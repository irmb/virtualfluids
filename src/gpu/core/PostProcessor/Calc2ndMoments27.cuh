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
//! \author Martin Schoenherr
//======================================================================================

#ifndef Calc2ndMoments27_CUH
#define Calc2ndMoments27_CUH

#include <basics/DataTypes.h>

void Calc2ndMomentsIncompSP27(real* kxyFromfcNEQ, real* kyzFromfcNEQ, real* kxzFromfcNEQ, real* kxxMyyFromfcNEQ,
                              real* kxxMzzFromfcNEQ, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                              unsigned int* neighborZ, unsigned long long numberOfLBnodes, unsigned int numberOfThreads,
                              real* DD, bool isEvenTimestep);

void Calc2ndMomentsCompSP27(real* kxyFromfcNEQ, real* kyzFromfcNEQ, real* kxzFromfcNEQ, real* kxxMyyFromfcNEQ,
                            real* kxxMzzFromfcNEQ, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                            unsigned int* neighborZ, unsigned long long numberOfLBnodes, unsigned int numberOfThreads,
                            real* DD, bool isEvenTimestep);

void Calc3rdMomentsIncompSP27(real* CUMbbb, real* CUMabc, real* CUMbac, real* CUMbca, real* CUMcba, real* CUMacb,
                              real* CUMcab, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                              unsigned int* neighborZ, unsigned long long numberOfLBnodes, unsigned int numberOfThreads,
                              real* DD, bool isEvenTimestep);

void Calc3rdMomentsCompSP27(real* CUMbbb, real* CUMabc, real* CUMbac, real* CUMbca, real* CUMcba, real* CUMacb, real* CUMcab,
                            unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                            unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD, bool isEvenTimestep);

void CalcHigherMomentsIncompSP27(real* CUMcbb, real* CUMbcb, real* CUMbbc, real* CUMcca, real* CUMcac, real* CUMacc,
                                 real* CUMbcc, real* CUMcbc, real* CUMccb, real* CUMccc, unsigned int* geoD,
                                 unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                                 unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD,
                                 bool isEvenTimestep);

void CalcHigherMomentsCompSP27(real* CUMcbb, real* CUMbcb, real* CUMbbc, real* CUMcca, real* CUMcac, real* CUMacc,
                               real* CUMbcc, real* CUMcbc, real* CUMccb, real* CUMccc, unsigned int* geoD,
                               unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                               unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD,
                               bool isEvenTimestep);

#endif

//! \}
