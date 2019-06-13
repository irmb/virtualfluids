#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#include <core/PointerDefinitions.h>

class Parameter;
class CudaMemoryManager;

void writeInit(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager);
void writeTimestep(Parameter* para, CudaMemoryManager* cudaManager, unsigned int t);
void writeParticle(Parameter* para, unsigned int t);

#endif
