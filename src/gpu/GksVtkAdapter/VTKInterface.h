#ifndef VTKInterface_H
#define VTKInterface_H

#include "GksGpu/Parameters/Parameters.h"

#include "VirtualFluidsDefinitions.h"

namespace GksGpu{ 
struct DataBase;
class TurbulenceAnalyzer;
struct ConcreteHeatFlux;
}

void VF_PUBLIC writeVtkXML(std::shared_ptr<GksGpu::DataBase> dataBase, 
                           GksGpu::Parameters parameters, 
                           int mode, 
                           std::string filename);

void VF_PUBLIC writeVtkXMLParallelSummaryFile(std::shared_ptr<GksGpu::DataBase> dataBase, 
                                              GksGpu::Parameters parameters, 
                                              std::string filename,
                                              uint mpiWorldSize);

void VF_PUBLIC writeTurbulenceVtkXML(std::shared_ptr<GksGpu::DataBase> dataBase, 
                                     std::shared_ptr<GksGpu::TurbulenceAnalyzer> turbulenceAnalyzer,
                                     int mode, 
                                     std::string filename);

void VF_PUBLIC writeTurbulenceVtkXMLParallelSummaryFile(std::shared_ptr<GksGpu::DataBase> dataBase, 
                                                        std::shared_ptr<GksGpu::TurbulenceAnalyzer> turbulenceAnalyzer,
                                                        GksGpu::Parameters parameters, 
                                                        std::string filename,
                                                        uint mpiWorldSize);

void VF_PUBLIC mapFlowField( std::shared_ptr<GksGpu::DataBase> base, std::shared_ptr<GksGpu::DataBase> target );

void VF_PUBLIC writeConcreteHeatFluxVtkXML(std::shared_ptr<GksGpu::DataBase> dataBase, 
                                           std::shared_ptr<GksGpu::ConcreteHeatFlux> bc, 
                                           GksGpu::Parameters parameters, 
                                           int mode, 
                                           std::string filename);

#endif