#ifndef CudaTimer_H
#define CudaTimer_H


#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>


class CudaTimer
{
public:
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
    double ftimeS = {0.0};

    cudaEvent_t stop_t;
    cudaEvent_t start_t;
    double ftimeE {0.0};

};

#endif
