#ifndef CalcTurbulenceIntensity_H
#define CalcTurbulenceIntensity_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint size);
extern "C" void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint size);
extern "C" void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint size);
extern "C" void resetTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint size);

#endif
