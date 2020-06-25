#ifndef DragLift_H
#define DragLift_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Utilities/StringUtil.hpp"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void calcDragLift(Parameter* para, CudaMemoryManager* cudaManager, int lev);
extern "C" void allocDragLift(Parameter* para, CudaMemoryManager* cudaManager);
extern "C" void printDragLift(Parameter* para, CudaMemoryManager* cudaManager, int timestep);

#endif
