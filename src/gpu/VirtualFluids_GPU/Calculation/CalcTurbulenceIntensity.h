#ifndef CalcTurbulenceIntensity_H
#define CalcTurbulenceIntensity_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint size);
extern "C" void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint size);
extern "C" void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint size);
extern "C" void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaManager, uint size);

void writeTurbulenceIntensityToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int size);
void writeVeloFluctuationToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int size);
void writeVeloMeansToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int size);
void writeAllTiDatafToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int size);

void writeTiStuffToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int size, std::vector<real *> &data,
                  std::vector<std::string> &datanames);

#endif
