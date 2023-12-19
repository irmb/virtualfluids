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
//! \addtogroup gpu_Communication Communication
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//======================================================================================

#ifndef EXCHANGEDATA27_DEVICE_CUH
#define EXCHANGEDATA27_DEVICE_CUH

#include <cuda.h>
#include <cuda_runtime.h>

#include <basics/DataTypes.h>

void GetSendFsPreDev27(real* DD, real* bufferFs, int* sendIndex, int buffmax, unsigned int* neighborX,
                       unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                       bool isEvenTimestep, unsigned int numberOfThreads, cudaStream_t stream = CU_STREAM_LEGACY);

void GetSendFsPostDev27(real* DD, real* bufferFs, int* sendIndex, int buffmax, unsigned int* neighborX,
                        unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                        bool isEvenTimestep, unsigned int numberOfThreads, cudaStream_t stream = CU_STREAM_LEGACY);

void SetRecvFsPreDev27(real* DD, real* bufferFs, int* recvIndex, int buffmax, unsigned int* neighborX,
                       unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                       bool isEvenTimestep, unsigned int numberOfThreads, cudaStream_t stream = CU_STREAM_LEGACY);

void SetRecvFsPostDev27(real* DD, real* bufferFs, int* recvIndex, int buffmax, unsigned int* neighborX,
                        unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                        bool isEvenTimestep, unsigned int numberOfThreads, cudaStream_t stream = CU_STREAM_LEGACY);

void getSendGsDevF3(real* G6, real* bufferGs, int* sendIndex, int buffmax, unsigned int* neighborX, unsigned int* neighborY,
                    unsigned int* neighborZ, unsigned long long numberOfLBnodes, bool isEvenTimestep,
                    unsigned int numberOfThreads);

void setRecvGsDevF3(real* G6, real* bufferGs, int* recvIndex, int buffmax, unsigned int* neighborX, unsigned int* neighborY,
                    unsigned int* neighborZ, unsigned long long numberOfLBnodes, bool isEvenTimestep,
                    unsigned int numberOfThreads);

#endif

//! \}
