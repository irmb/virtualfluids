#ifndef VTKInterface_H
#define VTKInterface_H

#include "GksGpu/DataBase/DataBase.h"
#include "GksGpu/Parameters/Parameters.h"

#include "VirtualFluidsDefinitions.h"

void VF_PUBLIC writeVtkXML(std::shared_ptr<DataBase> dataBase, 
                           Parameters parameters, 
                           int mode, 
                           std::string filename);


void VF_PUBLIC mapFlowField( std::shared_ptr<DataBase> base, std::shared_ptr<DataBase> target );

#endif