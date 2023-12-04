#ifndef CalcMean_H
#define CalcMean_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

void allocMean(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void allocMeanAD(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void calcMean(Parameter* para, unsigned int tdiff);
void resetMean(Parameter* para);

#endif
