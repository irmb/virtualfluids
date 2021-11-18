#ifndef FIND_TEMPERATURE_H
#define FIND_TEMPERATURE_H


class CudaMemoryManager;
class Parameter;


void initTemperatur(Parameter* para, CudaMemoryManager* cudaManager, int lev);

void findTempSim(Parameter* para, CudaMemoryManager* cudaManager);

void findTempVelSim(Parameter* para, CudaMemoryManager* cudaManager);

void findTempPressSim(Parameter* para, CudaMemoryManager* cudaManager);


#endif
