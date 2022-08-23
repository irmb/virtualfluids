#ifndef DEFINE_BCS_H
#define DEFINE_BCS_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

void findQ27(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void findBC27(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void findPressQShip(Parameter* para, CudaMemoryManager* cudaMemoryManager);

#endif
