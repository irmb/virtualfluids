#ifndef CalcTurbulenceIntensity_H
#define CalcTurbulenceIntensity_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint sizeOfTiArray);
extern "C" void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint sizeOfTiArray);
extern "C" void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint sizeOfTiArray);
extern "C" void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaManager, uint sizeOfTiArray);

void writeTurbulenceIntensityToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray);
void writeVeloFluctuationToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray);
void writeVeloMeansToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray);
void writeAllTiDatafToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray);

void writeTiStuffToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray, std::vector<real *> &data,
                  std::vector<std::string> &datanames);

#endif
