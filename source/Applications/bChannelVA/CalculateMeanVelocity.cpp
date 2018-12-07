#include "Postprocessing.h"
#include <vector>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTimerLog.h>
#include "MemoryUtil.h"

namespace PP
{
   void CalculateMeanVelocity(std::vector <std::string> dataName, std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_sum = vtkSmartPointer<vtkTimerLog>::New();
      timer_global->StartTimer();

      //mean velocity
      std::vector<double> Mvx;
      std::vector<double> Mvy;
      std::vector<double> Mvz;

      int numberOfPoints;
      vtkUnstructuredGrid* unstructuredGridOut = vtkUnstructuredGrid::New();

      int numberOfFilese = (int)dataName.size();

      print("calculate mean velocity: start");
      timer_calculate->StartTimer();

      for (int i = 0; i < numberOfFilese; i++)
      {
         timer_read->StartTimer();
         print("read data set from "+dataName[i]+": start");
         vtkDataSet* dataSet(ReadDataSet(dataName[i].c_str()));
         print("read data set from "+dataName[i]+": end");
         timer_read->StopTimer();
         print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

         vtkDataArray *vxArray = dataSet->GetPointData()->GetArray("Vx");
         vtkDataArray *vyArray = dataSet->GetPointData()->GetArray("Vy");
         vtkDataArray *vzArray = dataSet->GetPointData()->GetArray("Vz");

         if (i == 0)
         {
            numberOfPoints = dataSet->GetNumberOfPoints();
            Mvx.resize(numberOfPoints, 0.0);
            Mvy.resize(numberOfPoints, 0.0);
            Mvz.resize(numberOfPoints, 0.0);
            unstructuredGridOut->CopyStructure(dataSet);
         }

         timer_sum->StartTimer();
         print("calculate sum of velocities: start");

         for (int j = 0; j < numberOfPoints; j++)
         {
            Mvx[j] += vxArray->GetTuple1(j);
            Mvy[j] += vyArray->GetTuple1(j);
            Mvz[j] += vzArray->GetTuple1(j);
         }

         timer_sum->StopTimer();
         print("calculate sum of velocities: end");
         print("calculate sum of velocities time: "+toString(timer_sum->GetElapsedTime())+" s");

         dataSet->Delete();
         print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
      }

      for (int j = 0; j < numberOfPoints; j++)
      {
         Mvx[j] /= (double)numberOfFilese;
         Mvy[j] /= (double)numberOfFilese;
         Mvz[j] /= (double)numberOfFilese;
      }

      print("calculate mean velocity: end");
      timer_calculate->StopTimer();
      print("calculate mean velocity time: "+toString(timer_calculate->GetElapsedTime())+" s");

      std::string vtkfilename = output+".vtu";

      print("write data set to "+vtkfilename+": start");
      timer_write->StartTimer();

      vtkDoubleArray *mvxArray = vtkSmartPointer<vtkDoubleArray>::New();
      vtkDoubleArray *mvzArray = vtkSmartPointer<vtkDoubleArray>::New();
      vtkDoubleArray *mvyArray = vtkSmartPointer<vtkDoubleArray>::New();
      mvxArray->SetNumberOfComponents(1);
      mvyArray->SetNumberOfComponents(1);
      mvzArray->SetNumberOfComponents(1);

      mvxArray->SetName("Mvx");
      mvyArray->SetName("Mvy");
      mvzArray->SetName("Mvz");

      for (int i = 0; i < numberOfPoints; i++)
      {
         mvxArray->InsertNextValue(Mvx[i]);
         mvyArray->InsertNextValue(Mvy[i]);
         mvzArray->InsertNextValue(Mvz[i]);
      }

      unstructuredGridOut->GetPointData()->AddArray(mvxArray);
      unstructuredGridOut->GetPointData()->AddArray(mvyArray);
      unstructuredGridOut->GetPointData()->AddArray(mvzArray);


      vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
      writer->SetInputData(unstructuredGridOut);
      writer->SetFileName(vtkfilename.c_str());
      //writer->SetDataModeToBinary();
      //writer->SetDataModeToAscii();
      //writer->SetCompressorTypeToNone();
      writer->SetDataModeToAppended();
      writer->SetCompressorTypeToZLib();
      writer->Update();
      writer->Delete();

      mvxArray->Delete();
      mvyArray->Delete();
      mvzArray->Delete();

      unstructuredGridOut->Delete();

      print("write data set: end");
      timer_write->StopTimer();
      print("write data set time: "+toString(timer_write->GetElapsedTime())+" s");

      timer_global->StopTimer();
      print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
   }
}
