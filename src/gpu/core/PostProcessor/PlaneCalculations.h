#ifndef PLANE_CALCULATIONS_H
#define PLANE_CALCULATIONS_H

#include "Parameter/Parameter.h"
#include "Cuda/CudaMemoryManager.h"
#include "basics/utilities/UbSystem.h"

#include <iostream>
#include <stdio.h>

void setSizeOfPlane(Parameter* para, int lev, unsigned int z);
void calcPressure(Parameter* para, std::string inorout, int lev);
void calcFlowRate(Parameter* para, int lev);

//advection + diffusion
void calcPlaneConc(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
void allocPlaneConc(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void printPlaneConc(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void printRE(Parameter* para, CudaMemoryManager* cudaMemoryManager, int timestep);

#endif
