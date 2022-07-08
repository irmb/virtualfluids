#ifndef DEFINE_BCS_H
#define DEFINE_BCS_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void findQ27(Parameter* para, CudaMemoryManager* cudaMemoryManager);

extern "C" void findBC27(Parameter* para, CudaMemoryManager* cudaMemoryManager);

extern "C" void findPressQShip(Parameter* para, CudaMemoryManager* cudaMemoryManager);

#endif
