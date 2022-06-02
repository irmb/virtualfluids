
#include <iostream>
#include <cuda_runtime.h>
#include "UbScheduler.h"
#include "Timer.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"


void Timer::initTimer()
{
    cudaEventCreate(&this->start_t);
    cudaEventCreate(&this->stop_t );
}

void Timer::startTimer()
{ 
    checkCudaErrors(cudaEventRecord(this->start_t)); 
}

void Timer::stopTimer()
{
        checkCudaErrors(cudaEventRecord(this->stop_t));
        checkCudaErrors(cudaEventSynchronize(this->stop_t));
        checkCudaErrors(cudaEventElapsedTime(&this->elapsedTime, this->start_t, this->stop_t));
        this->totalElapsedTime += this->elapsedTime;
}

void Timer::resetTimer()
{
        this->elapsedTime = 0.0;
        this->totalElapsedTime = 0.0;
}

void Timer::outputPerformance(uint t, Parameter* para, vf::gpu::Communicator& communicator)
{
    real fnups      = 0.0;
    real bandwidth  = 0.0;
    
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        fnups       += 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP * pow(2.,lev) / (this->totalElapsedTime*1.0E6);
        bandwidth   += (27.0+1.0) * 4.0 * 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP  / (this->totalElapsedTime*1.0E9);
    }

    if(this->firstOutput && communicator.getPID() == 0) //only display the legend once
    {
        VF_LOG_INFO("PID \t --- {} ---  Processing time (ms) \t Nups in Mio \t Bandwidth in GB/sec", this->name );
        this->firstOutput = false;
    }

    VF_LOG_INFO(" {} \t --- {} --- {:>8.1f}/ {:<8.1f} \t   {:5.1f} \t       {:4.1f}",  communicator.getPID(), this->name, this->elapsedTime, this->totalElapsedTime, fnups, bandwidth);

    // When using multiple GPUs, sum the nups of all processes
    if (communicator.getNummberOfProcess() > 1) {
        double nupsSum =  communicator.sumNups(fnups);
        VF_LOG_INFO("Sum of all {}\t processes: Nups in Mio: {:.1f}", communicator.getNummberOfProcess(), nupsSum);
    }
}