#ifndef VTKAdapter_H
#define VTKAdapter_H

#include <vtkSmartPointer.h>
#include <vtkVersion.h>
 
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataObject.h>


#include <vtkIdList.h>

#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include <vtkXMLUnstructuredGridWriter.h>

#include <memory>
#include <functional>
#include <string>

#include "PointerDefinitions.h"

#include "VirtualFluidsDefinitions.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

typedef vtkSmartPointer<vtkUnstructuredGrid>          vtkGridPtr;
typedef vtkSmartPointer<vtkPoints>                    vtkPointsPtr;
typedef vtkSmartPointer<vtkIdList>                    vtkIdListPtr;
typedef vtkSmartPointer<vtkIntArray>                  vtkIntArrayPtr;
typedef vtkSmartPointer<vtkDoubleArray>               vtkDoubleArrayPtr;
typedef vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriterPtr;

struct rgbColor
{
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

vtkGridPtr VIRTUALFLUIDS_GPU_EXPORT getVtkUnstructuredOctGrid( SPtr<GksGpu::DataBase> dataBase, bool excludeGhostCells = false );

void VIRTUALFLUIDS_GPU_EXPORT addScalarIntCellData( vtkGridPtr grid,
                                     uint numberOfCells, 
                                     std::string name, 
                                     std::function<int(uint)> getData );

void VIRTUALFLUIDS_GPU_EXPORT addScalarRealCellData( vtkGridPtr grid,
                                      uint numberOfCells, 
                                      std::string name, 
                                      std::function<real(uint)> getData );

void VIRTUALFLUIDS_GPU_EXPORT addVectorCellData( vtkGridPtr grid,
                                  uint numberOfCells, 
                                  std::string name, 
                                  std::function<Vec3(uint)> getData );

void VIRTUALFLUIDS_GPU_EXPORT addBaseData( vtkGridPtr grid, SPtr<GksGpu::DataBase> dataBase, GksGpu::Parameters parameters );

void VIRTUALFLUIDS_GPU_EXPORT writeVtkUnstructuredGrid( vtkGridPtr grid, int mode, std::string filename );

void VIRTUALFLUIDS_GPU_EXPORT writeVtkParallelUnstructuredGridSummaryFile( vtkGridPtr grid, std::string filename, uint mpiWorldSize );

rgbColor VIRTUALFLUIDS_GPU_EXPORT colorMapCoolToWarmExtended( double value, double min, double max );

void VIRTUALFLUIDS_GPU_EXPORT writePNG( vtkDataObject* inputData, int nx, int ny, double L, double H, std::string filename );

#endif