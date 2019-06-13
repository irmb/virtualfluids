#ifndef VF_READER_H
#define VF_READER_H

#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "Input/kFullReader.h"
#include "Input/PositionReader.h"
//#include "GPU/GPU_Interface.h"

#include <iostream>

class CudaMemoryManager;

extern "C" void readVFkFull(Parameter* para, CudaMemoryManager* cudaManager, const std::string geometryFile);

extern "C" void readVFgeoFull(Parameter* para, const std::string geometryFile);

extern "C" void readVecSP(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readInterfaceCF(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readInterfaceFC(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readInterfaceOffCF(Parameter* para, CudaMemoryManager* cudaManager, const std::string geometryFile);

extern "C" void readInterfaceOffFC(Parameter* para, CudaMemoryManager* cudaManager, const std::string geometryFile);

extern "C" void readNoSlipBc(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readSlipBc(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readPressBc(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readPropellerCylinder(Parameter* para, CudaMemoryManager* cudaManager);

extern "C" void readMeasurePoints(Parameter* para, CudaMemoryManager* cudaManager);

#endif
