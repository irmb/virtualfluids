#ifndef PLANE_CALCULATIONS_H
#define PLANE_CALCULATIONS_H

#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"
#include "basics/utilities/UbSystem.h"

#include <iostream>
#include <stdio.h>

extern "C" void setSizeOfPlane(Parameter* para, int lev, unsigned int z);
extern "C" void calcPressure(Parameter* para, std::string inorout, int lev);
extern "C" void calcFlowRate(Parameter* para, int lev);

//advection + diffusion
extern "C" void calcPlaneConc(Parameter* para, CudaMemoryManager* cudaManager, int lev);
extern "C" void allocPlaneConc(Parameter* para, CudaMemoryManager* cudaManager);
extern "C" void printPlaneConc(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void printRE(Parameter* para, CudaMemoryManager* cudaManager, int timestep);

#endif
