#ifndef INIT_LATTICE_H
#define INIT_LATTICE_H

#include "Core/PointerDefinitions.h"

class Parameter;
class PreProcessor;
class CudaMemoryManager;

void initLattice(SPtr<Parameter> para, SPtr<PreProcessor> preProcessor, SPtr<CudaMemoryManager> cudaManager);

#endif
