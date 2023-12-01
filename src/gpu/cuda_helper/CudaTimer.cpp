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
