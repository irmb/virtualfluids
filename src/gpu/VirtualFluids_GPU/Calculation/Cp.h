#ifndef Cp_H
#define Cp_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

extern "C" void calcCp(Parameter* para, CudaMemoryManager* cudaManager, int lev);
extern "C" void printCpTopIntermediateStep(Parameter* para, unsigned int t, int lev);
extern "C" void printCpTop(Parameter* para, CudaMemoryManager* cudaManager, int lev);
extern "C" void printCpBottom(Parameter* para, CudaMemoryManager* cudaManager);
extern "C" void printCpBottom2(Parameter* para, CudaMemoryManager* cudaManager);



extern "C" void excludeGridInterfaceNodesForMirror(Parameter* para, int lev);
extern "C" void calcPressForMirror(Parameter* para, CudaMemoryManager* cudaManager, int lev);
//Ensight Gold
extern "C" void printCaseFile(Parameter* para);
extern "C" void printGeoFile(Parameter* para, bool fileFormat);
extern "C" void printScalars(Parameter* para, bool fileFormat);
//functions to write binary files
extern "C" void writeIntToFile(const int &i, std::ofstream &ofile);
extern "C" void writeFloatToFile(const float &f, std::ofstream &ofile);
extern "C" void writeStringToFile(const std::string &s, std::ofstream &ofile);

#endif
