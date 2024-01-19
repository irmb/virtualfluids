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
//! \addtogroup gpu_cuda_helper cuda_helper
//! \{
//! \author Soeren Peters
//=======================================================================================
#include "CudaTimer.h"

namespace vf::cuda
{

void CudaTimer::createSdkTimer()
{
    sdkCreateTimer(&sdkTimer);
}

void CudaTimer::startSdkTimer()
{
    sdkStartTimer(&sdkTimer);
}

void CudaTimer::createEventTimer()
{
    checkCudaErrors(cudaEventCreate(&start_t));
    checkCudaErrors(cudaEventCreate(&stop_t));
}

void CudaTimer::startEventTimer()
{
    checkCudaErrors(cudaEventRecord(start_t));
}

void CudaTimer::stopSdkTimer(float &timeSinceLastStop,double &totalTime)
{
    sdkStopTimer(&sdkTimer);
    timeSinceLastStop = sdkGetTimerValue(&sdkTimer);
    sdkResetTimer(&sdkTimer);
    ftimeS += static_cast<double>(timeSinceLastStop);
    totalTime = ftimeS;
}

void CudaTimer::stopEventTimer(float &timeSinceLastStop,double &totalTime)
{
    checkCudaErrors(cudaEventRecord(stop_t));
    checkCudaErrors(cudaEventSynchronize(stop_t));
    checkCudaErrors(cudaEventElapsedTime(&timeSinceLastStop, start_t, stop_t));
    ftimeE += static_cast<double>(timeSinceLastStop);
    totalTime = ftimeE;
}

void CudaTimer::deleteSdkTimer()
{
    sdkDeleteTimer(&sdkTimer);
}

void CudaTimer::deleteEventTimer()
{
    checkCudaErrors(cudaEventDestroy(start_t));
    checkCudaErrors(cudaEventDestroy(stop_t));
}

} // namespace vf::cuda

//! \}
