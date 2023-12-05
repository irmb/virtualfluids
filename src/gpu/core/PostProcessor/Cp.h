#ifndef Cp_H
#define Cp_H

#include "Calculation/Calculation.h"

#include "Parameter/Parameter.h"
#include "Cuda/CudaMemoryManager.h"

void calcCp(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
void printCpTopIntermediateStep(Parameter* para, unsigned int t, int lev);
void printCpTop(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
void printCpBottom(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void printCpBottom2(Parameter* para, CudaMemoryManager* cudaMemoryManager);



void excludeGridInterfaceNodesForMirror(Parameter* para, int lev);
void calcPressForMirror(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
//Ensight Gold
void printCaseFile(Parameter* para);
void printGeoFile(Parameter* para, bool fileFormat);
void printScalars(Parameter* para, bool fileFormat);
//functions to write binary files
void writeIntToFile(const int &i, std::ofstream &ofile);
void writeFloatToFile(const float &f, std::ofstream &ofile);
void writeStringToFile(const std::string &s, std::ofstream &ofile);

#endif
