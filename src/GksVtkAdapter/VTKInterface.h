#ifndef VTKInterface_H
#define VTKInterface_H

#include "GksGpu/Parameters/Parameters.h"

#include "VirtualFluidsDefinitions.h"

struct DataBase;
class TurbulenceAnalyzer;

void VF_PUBLIC writeVtkXML(std::shared_ptr<DataBase> dataBase, 
                           Parameters parameters, 
                           int mode, 
                           std::string filename);

void VF_PUBLIC writeTurbulenceVtkXML(std::shared_ptr<DataBase> dataBase, 
                                     std::shared_ptr<TurbulenceAnalyzer> turbulenceAnalyzer,
                                     int mode, 
                                     std::string filename);

void VF_PUBLIC mapFlowField( std::shared_ptr<DataBase> base, std::shared_ptr<DataBase> target );

#endif