#ifndef TIMER_H
#define TIMER_H
#include <cuda_runtime.h>

#include "DataTypes.h"
#include "Parameter/Parameter.h"
#include <logger/Logger.h>

namespace vf::parallel
{
class Communicator;
}
class Parameter;

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
    void outputPerformance(uint t, Parameter* para, vf::parallel::Communicator& communicator);

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