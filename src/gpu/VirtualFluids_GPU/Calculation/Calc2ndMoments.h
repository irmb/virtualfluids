#ifndef Calc2ndMoments_H
#define Calc2ndMoments_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

//2nd
void alloc2ndMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void init2ndMoments(Parameter* para);
void calc2ndMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);

//3rd
void alloc3rdMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void init3rdMoments(Parameter* para);
void calc3rdMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);

//higher order
void allocHigherOrderMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void initHigherOrderMoments(Parameter* para);
void calcHigherOrderMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);

#endif
