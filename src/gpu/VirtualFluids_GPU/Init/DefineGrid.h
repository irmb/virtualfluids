#ifndef DEFINE_GRID_H
#define DEFINE_GRID_H

#include "LBM/LB.h"
#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void defineGrid(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager);

#endif
