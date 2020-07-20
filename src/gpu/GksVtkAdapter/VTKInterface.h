#ifndef VTKInterface_H
#define VTKInterface_H

#include "GksGpu/Parameters/Parameters.h"

#include "VirtualFluidsDefinitions.h"

namespace GksGpu{ 
struct DataBase;
class TurbulenceAnalyzer;
struct ConcreteHeatFlux;
}

void VIRTUALFLUIDS_GPU_EXPORT writeVtkXML(std::shared_ptr<GksGpu::DataBase> dataBase,
                           GksGpu::Parameters parameters, 
                           int mode, 
                           std::string filename);

void VIRTUALFLUIDS_GPU_EXPORT writeVtkXMLParallelSummaryFile(std::shared_ptr<GksGpu::DataBase> dataBase,
                                              GksGpu::Parameters parameters, 
                                              std::string filename,
                                              uint mpiWorldSize);

void VIRTUALFLUIDS_GPU_EXPORT writeTurbulenceVtkXML(std::shared_ptr<GksGpu::DataBase> dataBase,
                                     std::shared_ptr<GksGpu::TurbulenceAnalyzer> turbulenceAnalyzer,
                                     int mode, 
                                     std::string filename);

void VIRTUALFLUIDS_GPU_EXPORT writeTurbulenceVtkXMLParallelSummaryFile(std::shared_ptr<GksGpu::DataBase> dataBase,
                                                        std::shared_ptr<GksGpu::TurbulenceAnalyzer> turbulenceAnalyzer,
                                                        GksGpu::Parameters parameters, 
                                                        std::string filename,
                                                        uint mpiWorldSize);

void VIRTUALFLUIDS_GPU_EXPORT mapFlowField( std::shared_ptr<GksGpu::DataBase> base, std::shared_ptr<GksGpu::DataBase> target );

void VIRTUALFLUIDS_GPU_EXPORT writeConcreteHeatFluxVtkXML(std::shared_ptr<GksGpu::DataBase> dataBase,
                                           std::shared_ptr<GksGpu::ConcreteHeatFlux> bc, 
                                           GksGpu::Parameters parameters, 
                                           int mode, 
                                           std::string filename);

#endif