#ifndef Particles_H
#define Particles_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "GPU/CudaMemoryManager.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Parameter/Parameter.h"

//void calcDragLift(Parameter* para, int lev);
void allocParticles(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void initParticles(Parameter* para);
void propagateParticles(Parameter* para, unsigned int t);
void copyAndPrintParticles(Parameter* para, CudaMemoryManager* cudaMemoryManager, unsigned int t, bool isInit);

void rearrangeGeometry(Parameter* para, CudaMemoryManager* cudaMemoryManager);

#endif
