#ifndef TIMER_H
#define TIMER_H

#include "helper_cuda.h"
#include <cuda_runtime.h>
#include "Core/DataTypes.h"

#include "UbScheduler.h"
#include "logger/Logger.h"
#include "Parameter/Parameter.h"

namespace vf::gpu{
    class Communicator;
}

class Timer
{
    public:
    Timer(std::string _name): name(_name)
    {
        this->initTimer();
    };
    
    ~Timer()
    {
        cudaEventDestroy(this->start_t);
        cudaEventDestroy(this->stop_t);
    };

    void initTimer();
    void startTimer();
    void stopTimer();
    void resetTimer();
    void outputPerformance(uint t, Parameter* para, vf::gpu::Communicator& communicator);

    float getElapsedTime(){ return this->elapsedTime; }
    float getTotalElapsedTime(){ return this->totalElapsedTime; }

    private:
    
    cudaEvent_t start_t, stop_t;
    float elapsedTime = 0.0;
    float totalElapsedTime = 0.0;
    std::string name;

    bool firstOutput = true;
};



#endif 