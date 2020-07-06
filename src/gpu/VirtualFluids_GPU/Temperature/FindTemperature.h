#ifndef FIND_TEMPERATURE_H
#define FIND_TEMPERATURE_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Temperature/FindQTemp.h"

class CudaMemoryManager;

extern "C" void initTemperatur(Parameter* para, CudaMemoryManager* cudaManager, int lev);

extern "C" void findTempSim(Parameter* para, CudaMemoryManager* cudaManager);
							
extern "C" void findTempVelSim(Parameter* para, CudaMemoryManager* cudaManager);
								
extern "C" void findTempPressSim(Parameter* para, CudaMemoryManager* cudaManager);

#endif
