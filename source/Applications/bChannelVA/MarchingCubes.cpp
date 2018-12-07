#include <vtkImageReader.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkSTLWriter.h>
#include <vtkVectorNorm.h>
#include <vtkFileOutputWindow.h>
#include <vtkMergePoints.h>

void MarchingCubes()
{
   vtkSmartPointer<vtkFileOutputWindow> fileOutputWindow = vtkSmartPointer<vtkFileOutputWindow>::New();
   fileOutputWindow->SetFileName("D:/Projects/SFB880/GeometrienPoroeseMedien/Membran/ow.txt");
   fileOutputWindow->SetFlush(0);
   fileOutputWindow->SetInstance(fileOutputWindow);

   vtkSmartPointer<vtkMergePoints> locator = vtkSmartPointer<vtkMergePoints>::New();
   locator->SetDivisions(184.5, 178, 50);
   locator->SetNumberOfPointsPerBucket(2);
   locator->AutomaticOff();

   std::string filePath = "D:/Projects/SFB880/GeometrienPoroeseMedien/Membran/membran370x357x101.raw";
   // Particles.raw supplied by VTK is big endian encoded
   //std::string filePath = "C://VTK//vtkdata-5.8.0//Data//Particles.raw";
   // Read the file
   vtkSmartPointer<vtkImageReader> reader = vtkSmartPointer<vtkImageReader>::New();
   reader->SetFileName(filePath.c_str());
   reader->SetDataByteOrderToBigEndian();
   reader->SetDataScalarTypeToChar();
   reader->SetDataExtent(0, 369, 0, 356, 0, 100);
   reader->SetDataSpacing(1, 1, 1);
   reader->Update();

   vtkSmartPointer<vtkMarchingCubes> iso = vtkSmartPointer<vtkMarchingCubes>::New();
   iso->SetInputConnection(reader->GetOutputPort());
   iso->SetValue(0, 62);
   iso->ComputeNormalsOff();
   iso->ComputeScalarsOn();
   iso->ComputeGradientsOff();
   iso->SetLocator(locator);

   //vtkSmartPointer<vtkVectorNorm> gradient = vtkSmartPointer<vtkVectorNorm>::New();
   //gradient->SetInputConnection(iso->GetOutputPort());

   vtkSmartPointer<vtkSTLWriter> stl = vtkSmartPointer<vtkSTLWriter>::New();
   stl->SetInputConnection(iso->GetOutputPort());
   stl->SetFileName("D:/Projects/SFB880/GeometrienPoroeseMedien/Membran/membran.stl");
   stl->Write();

}