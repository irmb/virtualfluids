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
   void CalculateFluctuations(std::vector <std::string> dataNames, std::string meanVelocityData, std::string output)
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
      std::vector<double> Vxy;
      std::vector<double> Vxz;
      std::vector<double> Vyz;

      std::vector<double> Vxxx;
      std::vector<double> Vyyy;
      std::vector<double> Vzzz;
      std::vector<double> Vxxy;
      std::vector<double> Vxxz;
      std::vector<double> Vyyx;
      std::vector<double> Vyyz;
      std::vector<double> Vzzx;
      std::vector<double> Vzzy;
      std::vector<double> Vxyz;

      int numberOfPoints;
      vtkUnstructuredGrid *unstructuredGridOut = vtkUnstructuredGrid::New();

      int numberOfFilese = (int)dataNames.size();

      print("read data set from "+meanVelocityData+": start");
      timer_read->StartTimer();
      vtkDataSet *dataSetMv(ReadDataSet(meanVelocityData.c_str()));
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

         if (i == 0)
         {
            numberOfPoints = dataSet->GetNumberOfPoints();
            Vxx.resize(numberOfPoints, 0.0);
            Vyy.resize(numberOfPoints, 0.0);
            Vzz.resize(numberOfPoints, 0.0);
            Vxy.resize(numberOfPoints, 0.0);
            Vxz.resize(numberOfPoints, 0.0);
            Vyz.resize(numberOfPoints, 0.0);

            Vxxx.resize(numberOfPoints, 0.0);
            Vyyy.resize(numberOfPoints, 0.0);
            Vzzz.resize(numberOfPoints, 0.0);
            Vxxy.resize(numberOfPoints, 0.0);
            Vxxz.resize(numberOfPoints, 0.0);
            Vyyx.resize(numberOfPoints, 0.0);
            Vyyz.resize(numberOfPoints, 0.0);
            Vzzx.resize(numberOfPoints, 0.0);
            Vzzy.resize(numberOfPoints, 0.0);
            Vxyz.resize(numberOfPoints, 0.0);
            unstructuredGridOut->CopyStructure(dataSet);
         }

         timer_sum->StartTimer();
         print("calculate sum of velocities: start");
         for (int j = 0; j < numberOfPoints; j++)
         {
            double tfx = vxArray->GetTuple1(j) - mvxArray->GetTuple1(j);
            double tfy = vyArray->GetTuple1(j) - mvyArray->GetTuple1(j);
            double tfz = vzArray->GetTuple1(j) - mvzArray->GetTuple1(j);

            Vxx[j] += tfx*tfx;
            Vyy[j] += tfy*tfy;
            Vzz[j] += tfz*tfz;
            Vxy[j] += tfx*tfy;
            Vxz[j] += tfx*tfz;
            Vyz[j] += tfy*tfz;

            Vxxx[j] += tfx*tfx*tfx;
            Vyyy[j] += tfy*tfy*tfy;
            Vzzz[j] += tfz*tfz*tfz;
            Vxxy[j] += tfx*tfx*tfy;
            Vxxz[j] += tfx*tfx*tfz;
            Vyyx[j] += tfy*tfy*tfx;
            Vyyz[j] += tfy*tfy*tfz;
            Vzzx[j] += tfz*tfz*tfx;
            Vzzy[j] += tfz*tfz*tfy;
            Vxyz[j] += tfx*tfy*tfz;
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
         Vxy[j] = Vxy[j]/(double)numberOfFilese;
         Vxz[j] = Vxz[j]/(double)numberOfFilese;
         Vyz[j] = Vyz[j]/(double)numberOfFilese;

         Vxxx[j] = Vxxx[j]/(double)numberOfFilese;
         Vyyy[j] = Vyyy[j]/(double)numberOfFilese;
         Vzzz[j] = Vzzz[j]/(double)numberOfFilese;
         Vxxy[j] = Vxxy[j]/(double)numberOfFilese;
         Vxxz[j] = Vxxz[j]/(double)numberOfFilese;
         Vyyx[j] = Vyyx[j]/(double)numberOfFilese;
         Vyyz[j] = Vyyz[j]/(double)numberOfFilese;
         Vzzx[j] = Vzzx[j]/(double)numberOfFilese;
         Vzzy[j] = Vzzy[j]/(double)numberOfFilese;
         Vxyz[j] = Vxyz[j]/(double)numberOfFilese;
      }
   }

   print("calculate rms: end");
   timer_calculate->StopTimer();
   print("calculate rms time: "+toString(timer_calculate->GetElapsedTime())+" s");

   std::string vtkfilename = output+".vtu";

   print("write data set to "+vtkfilename+": start");
   timer_write->StartTimer();

   vtkDoubleArray *vxxArray = vtkDoubleArray::New();
   vtkDoubleArray *vyyArray = vtkDoubleArray::New();
   vtkDoubleArray *vzzArray = vtkDoubleArray::New();
   vtkDoubleArray *vxyArray = vtkDoubleArray::New();
   vtkDoubleArray *vxzArray = vtkDoubleArray::New();
   vtkDoubleArray *vyzArray = vtkDoubleArray::New();
   vtkDoubleArray *vxxxArray = vtkDoubleArray::New();
   vtkDoubleArray *vyyyArray = vtkDoubleArray::New();
   vtkDoubleArray *vzzzArray = vtkDoubleArray::New();
   vtkDoubleArray *vxxyArray = vtkDoubleArray::New();
   vtkDoubleArray *vxxzArray = vtkDoubleArray::New();
   vtkDoubleArray *vyyxArray = vtkDoubleArray::New();
   vtkDoubleArray *vyyzArray = vtkDoubleArray::New();
   vtkDoubleArray *vzzxArray = vtkDoubleArray::New();
   vtkDoubleArray *vzzyArray = vtkDoubleArray::New();
   vtkDoubleArray *vxyzArray = vtkDoubleArray::New();

   vxxArray->SetNumberOfComponents(1);
   vyyArray->SetNumberOfComponents(1);
   vzzArray->SetNumberOfComponents(1);
   vxyArray->SetNumberOfComponents(1);
   vxzArray->SetNumberOfComponents(1);
   vyzArray->SetNumberOfComponents(1);

   vxxxArray->SetNumberOfComponents(1);
   vyyyArray->SetNumberOfComponents(1);
   vzzzArray->SetNumberOfComponents(1);
   vxxyArray->SetNumberOfComponents(1);
   vxxzArray->SetNumberOfComponents(1);
   vyyxArray->SetNumberOfComponents(1);
   vyyzArray->SetNumberOfComponents(1);
   vzzxArray->SetNumberOfComponents(1);
   vzzyArray->SetNumberOfComponents(1);
   vxyzArray->SetNumberOfComponents(1);

   vxxArray->SetName("Vxx");
   vyyArray->SetName("Vyy");
   vzzArray->SetName("Vzz");
   vxyArray->SetName("Vxy");
   vxzArray->SetName("Vxz");
   vyzArray->SetName("Vyz");

   vxxxArray->SetName("Vxxx");
   vyyyArray->SetName("Vyyy");
   vzzzArray->SetName("Vzzz");
   vxxyArray->SetName("Vxxy");
   vxxzArray->SetName("Vxxz");
   vyyxArray->SetName("Vyyx");
   vyyzArray->SetName("Vyyz");
   vzzxArray->SetName("Vzzx");
   vzzyArray->SetName("Vzzy");
   vxyzArray->SetName("Vxyz");


   for (int i = 0; i < numberOfPoints; i++)
   {
      vxxArray->InsertNextValue(Vxx[i]);
      vyyArray->InsertNextValue(Vyy[i]);
      vzzArray->InsertNextValue(Vzz[i]);
      vxyArray->InsertNextValue(Vxy[i]);
      vxzArray->InsertNextValue(Vxz[i]);
      vyzArray->InsertNextValue(Vyz[i]);

      vxxxArray->InsertNextValue(Vxxx[i]);
      vyyyArray->InsertNextValue(Vyyy[i]);
      vzzzArray->InsertNextValue(Vzzz[i]);
      vxxyArray->InsertNextValue(Vxxy[i]);
      vxxzArray->InsertNextValue(Vxxz[i]);
      vyyxArray->InsertNextValue(Vyyx[i]);
      vyyzArray->InsertNextValue(Vyyz[i]);
      vzzxArray->InsertNextValue(Vzzx[i]);
      vzzyArray->InsertNextValue(Vzzy[i]);
      vxyzArray->InsertNextValue(Vxyz[i]);
   }

   unstructuredGridOut->GetPointData()->AddArray(vxxArray);
   unstructuredGridOut->GetPointData()->AddArray(vyyArray);
   unstructuredGridOut->GetPointData()->AddArray(vzzArray);
   unstructuredGridOut->GetPointData()->AddArray(vxyArray);
   unstructuredGridOut->GetPointData()->AddArray(vxzArray);
   unstructuredGridOut->GetPointData()->AddArray(vyzArray);

   unstructuredGridOut->GetPointData()->AddArray(vxxxArray);
   unstructuredGridOut->GetPointData()->AddArray(vyyyArray);
   unstructuredGridOut->GetPointData()->AddArray(vzzzArray);
   unstructuredGridOut->GetPointData()->AddArray(vxxyArray);
   unstructuredGridOut->GetPointData()->AddArray(vxxzArray);
   unstructuredGridOut->GetPointData()->AddArray(vyyxArray);
   unstructuredGridOut->GetPointData()->AddArray(vyyzArray);
   unstructuredGridOut->GetPointData()->AddArray(vzzxArray);
   unstructuredGridOut->GetPointData()->AddArray(vzzyArray);
   unstructuredGridOut->GetPointData()->AddArray(vxyzArray);

   vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
   writer->SetInputData(unstructuredGridOut);
   writer->SetFileName(vtkfilename.c_str());
   writer->SetDataModeToAppended();
   writer->SetCompressorTypeToZLib();
   writer->Update();
   writer->Delete();

   unstructuredGridOut->Delete();

   vxxArray ->Delete();
   vyyArray ->Delete();
   vzzArray ->Delete();
   vxyArray ->Delete();
   vxzArray ->Delete();
   vyzArray ->Delete();
   vxxxArray->Delete();
   vyyyArray->Delete();
   vzzzArray->Delete();
   vxxyArray->Delete();
   vxxzArray->Delete();
   vyyxArray->Delete();
   vyyzArray->Delete();
   vzzxArray->Delete();
   vzzyArray->Delete();
   vxyzArray->Delete();

   dataSetMv->Delete();

   print("write data set: end");
   timer_write->StopTimer();
   print("write data set time: "+toString(timer_write->GetElapsedTime())+" s");

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
   }
}
