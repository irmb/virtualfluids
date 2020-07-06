#ifndef CudaTimer_H
#define CudaTimer_H


#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>


class CudaTimer
{
public:
    CudaTimer();

    void createSdkTimer();
    void startSdkTimer();
    void stopSdkTimer(float &timeSinceLastStop,double &totalTime);
    void deleteSdkTimer();

    void createEventTimer();
    void startEventTimer();
    void stopEventTimer(float &timeSinceLastStop,double &totalTime);
    void deleteEventTimer();

private:
    StopWatchInterface *sdkTimer;
    double ftimeS;

    cudaEvent_t stop_t;
    cudaEvent_t start_t;
    double ftimeE;

};

#endif
