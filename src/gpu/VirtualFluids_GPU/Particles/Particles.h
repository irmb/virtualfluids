#ifndef Particles_H
#define Particles_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "GPU/CudaMemoryManager.h"
#include "Utilities/StringUtil.hpp"
#include "Parameter/Parameter.h"

//extern "C" void calcDragLift(Parameter* para, int lev);
extern "C" void allocParticles(Parameter* para, CudaMemoryManager* cudaManager);
extern "C" void initParticles(Parameter* para);
extern "C" void propagateParticles(Parameter* para, unsigned int t);
extern "C" void copyAndPrintParticles(Parameter* para, CudaMemoryManager* cudaManager, unsigned int t, bool isInit);

extern "C" void rearrangeGeometry(Parameter* para, CudaMemoryManager* cudaManager);

#endif