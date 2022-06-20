#ifndef VF_READER_H
#define VF_READER_H


class CudaMemoryManager;
class Parameter;


void readPropellerCylinder(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void readMeasurePoints(Parameter* para, CudaMemoryManager* cudaMemoryManager);

#endif
