#ifndef DragLift_H
#define DragLift_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

void calcDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
void allocDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void printDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager, int timestep);

#endif
