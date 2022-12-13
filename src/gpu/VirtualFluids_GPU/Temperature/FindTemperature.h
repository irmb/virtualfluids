#ifndef FIND_TEMPERATURE_H
#define FIND_TEMPERATURE_H


class CudaMemoryManager;
class Parameter;


void initTemperatur(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);

void findTempSim(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void findTempVelSim(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void findTempPressSim(Parameter* para, CudaMemoryManager* cudaMemoryManager);


#endif
