#ifndef Calc2ndMoments_H
#define Calc2ndMoments_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

//2nd
extern "C" void alloc2ndMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);
extern "C" void init2ndMoments(Parameter* para);
extern "C" void calc2ndMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);

//3rd
extern "C" void alloc3rdMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);
extern "C" void init3rdMoments(Parameter* para);
extern "C" void calc3rdMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);

//higher order
extern "C" void allocHigherOrderMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);
extern "C" void initHigherOrderMoments(Parameter* para);
extern "C" void calcHigherOrderMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager);

#endif
