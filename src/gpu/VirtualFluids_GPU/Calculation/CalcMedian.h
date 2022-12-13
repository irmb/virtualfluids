#ifndef CalcMedian_H
#define CalcMedian_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

void allocMedian(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void allocMedianAD(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void calcMedian(Parameter* para, unsigned int tdiff);
void resetMedian(Parameter* para);

#endif
