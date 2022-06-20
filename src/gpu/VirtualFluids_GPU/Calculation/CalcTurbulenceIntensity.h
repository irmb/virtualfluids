#ifndef CalcTurbulenceIntensity_H
#define CalcTurbulenceIntensity_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaMemoryManager);
extern "C" void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaMemoryManager, uint tdiff);
extern "C" void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaMemoryManager, uint tdiff);
extern "C" void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaMemoryManager);
extern "C" void cudaFreeTurbulenceIntensityArrays(Parameter *para, CudaMemoryManager *cudaMemoryManager);


void writeTurbulenceIntensityToFile(Parameter *para, uint timestep);
void writeVeloFluctuationToFile(Parameter *para, uint timeste);
void writeVeloMeansToFile(Parameter *para, uint timestep);
void writeAllTiDatafToFile(Parameter *para, uint timestep);

void writeTiStuffToFile(Parameter *para, uint timestep, int sizeOfTiArray, std::vector<real *> &data,
                  std::vector<std::string> &datanames);

#endif
