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

void CalculateRMS(std::vector <std::string> dataNames, std::string meanVelocityData, std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_sum = vtkSmartPointer<vtkTimerLog>::New();
   timer_global->StartTimer();

   //mean velocity
   std::vector<double> Vxx;
   std::vector<double> Vyy;
   std::vector<double> Vzz;

   int numberOfPoints;
   vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridOut = vtkSmartPointer<vtkUnstructuredGrid>::New();

   int numberOfFilese = (int)dataNames.size();

   print("read data set from "+meanVelocityData+": start");
   timer_read->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSetMv(ReadDataSet(meanVelocityData.c_str()));
   print("read data set from "+meanVelocityData+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   vtkSmartPointer<vtkDataArray> mvxArray = dataSetMv->GetPointData()->GetArray("Mvx");
   vtkSmartPointer<vtkDataArray> mvyArray = dataSetMv->GetPointData()->GetArray("Mvy");
   vtkSmartPointer<vtkDataArray> mvzArray = dataSetMv->GetPointData()->GetArray("Mvz");

   print("calculate mean velocity: start");
   timer_calculate->StartTimer();

   for (int i = 0; i < numberOfFilese; i++)
   {
      timer_read->StartTimer();
      print("read data set from "+dataNames[i]+": start");
      vtkDataSet* dataSet(ReadDataSet(dataNames[i].c_str()));
      print("read data set from "+dataNames[i]+": end");
      timer_read->StopTimer();
      print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

      vtkSmartPointer<vtkDataArray> vxArray = dataSet->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSet->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSet->GetPointData()->GetArray("Vz");

      vtkDoubleArray* vxxArray = vtkDoubleArray::New();
      vxxArray->SetNumberOfComponents(1);
      vxxArray->SetName( "Vxx" );

      vtkDoubleArray* vyyArray = vtkDoubleArray::New();
      vyyArray->SetNumberOfComponents(1);
      vyyArray->SetName( "Vyy" );

      vtkDoubleArray* vzzArray = vtkDoubleArray::New();
      vzzArray->SetNumberOfComponents(1);
      vzzArray->SetName( "Vzz" );

      if (i == 0)
      {
         numberOfPoints = dataSet->GetNumberOfPoints();
         Vxx.resize(numberOfPoints,0.0);
         Vyy.resize(numberOfPoints,0.0);
         Vzz.resize(numberOfPoints,0.0);
         unstructuredGridOut->CopyStructure(dataSet);
      }

      timer_sum->StartTimer();
      print("calculate sum of velocities: start");
      for (int j = 0; j < numberOfPoints; j++)
      {
         double vxx = vxArray->GetTuple1(j) - mvxArray->GetTuple1(j);
         double vyy = vyArray->GetTuple1(j) - mvyArray->GetTuple1(j);
         double vzz = vzArray->GetTuple1(j) - mvzArray->GetTuple1(j);
         Vxx[j] += sqrt(vxx*vxx);
         Vyy[j] += sqrt(vyy*vyy);
         Vzz[j] += sqrt(vzz*vzz);
         vxxArray->InsertNextValue(Vxx[j]);
         vyyArray->InsertNextValue(Vyy[j]);
         vzzArray->InsertNextValue(Vzz[j]);
      }

      timer_sum->StopTimer();
      print("calculate sum of velocities: end");
      print("calculate sum of velocities time: "+toString(timer_sum->GetElapsedTime())+" s");

      dataSet->Delete();
      print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe()/1e9)+" GByte");
   }

   #pragma omp parallel num_threads(8)
   {
      #pragma omp for 
      for (int j = 0; j < numberOfPoints; j++)
      {
         Vxx[j] = Vxx[j]/(double)numberOfFilese;
         Vyy[j] = Vyy[j]/(double)numberOfFilese;
         Vzz[j] = Vzz[j]/(double)numberOfFilese;
      }
   }

   print("calculate rms: end");
   timer_calculate->StopTimer();
   print("calculate rms time: "+toString(timer_calculate->GetElapsedTime())+" s");

   std::string vtkfilename = output+".vtu";

   print("write data set to "+vtkfilename+": start");
   timer_write->StartTimer();

   vtkSmartPointer<vtkDoubleArray> vxxArray = vtkSmartPointer<vtkDoubleArray>::New();
   vxxArray->SetNumberOfComponents(1);
   vxxArray->SetName( "Vxx" );

   vtkSmartPointer<vtkDoubleArray> vyyArray = vtkSmartPointer<vtkDoubleArray>::New();
   vyyArray->SetNumberOfComponents(1);
   vyyArray->SetName( "Vyy" );

   vtkSmartPointer<vtkDoubleArray> vzzArray = vtkSmartPointer<vtkDoubleArray>::New();
   vzzArray->SetNumberOfComponents(1);
   vzzArray->SetName( "Vzz" );

   for(int i = 0; i < numberOfPoints; i++)
   {
      vxxArray->InsertNextValue(Vxx[i]);
      vyyArray->InsertNextValue(Vyy[i]);
      vzzArray->InsertNextValue(Vzz[i]);
   }

   unstructuredGridOut->GetPointData()->AddArray(vxxArray);
   unstructuredGridOut->GetPointData()->AddArray(vyyArray);
   unstructuredGridOut->GetPointData()->AddArray(vzzArray);

   vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
   writer->SetInputData(unstructuredGridOut);
   writer->SetFileName(vtkfilename.c_str());
   writer->SetDataModeToAscii();
   writer->SetCompressorTypeToNone();
   writer->Update();

   print("write data set: end");
   timer_write->StopTimer();
   print("write data set time: "+toString(timer_write->GetElapsedTime())+" s");

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}

