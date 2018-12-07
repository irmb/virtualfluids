#include "Postprocessing.h"
#include <vector>
#include <vtkDataArray.h>
#include <vtkKdTreePointLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkTimerLog.h>
#include "MemoryUtil.h"

void VolumeAveraging(std::string dataName, double radius, std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_locator = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();
   timer_global->StartTimer();

   //velocity
   std::vector<double> Mvx;
   std::vector<double> Mvy;
   std::vector<double> Mvz;
   
   print("read data set from "+dataName+": start");
   timer_read->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSet(ReadDataSet(dataName.c_str()));
   int numberOfPoints = dataSet->GetNumberOfPoints();
   print("read data set from "+dataName+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   Mvx.resize(numberOfPoints,0.0);
   Mvy.resize(numberOfPoints,0.0);
   Mvz.resize(numberOfPoints,0.0);

   vtkSmartPointer<vtkDataArray> vxArray = dataSet->GetPointData()->GetArray("Mvx");
   vtkSmartPointer<vtkDataArray> vyArray = dataSet->GetPointData()->GetArray("Mvy");
   vtkSmartPointer<vtkDataArray> vzArray = dataSet->GetPointData()->GetArray("Mvz");

   print("build locator: start");
   timer_locator->StartTimer();
   // Create the tree
   vtkSmartPointer<vtkKdTreePointLocator> pointLocator = vtkSmartPointer<vtkKdTreePointLocator>::New();
   pointLocator->SetDataSet(dataSet);
   pointLocator->AutomaticOn();
   //pointLocator->SetNumberOfPointsPerBucket(2000);
   pointLocator->BuildLocator();
   timer_locator->StopTimer();
   print("build locator: end");
   print("build locator time: "+toString(timer_locator->GetElapsedTime())+" s");

   vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
   double center[3];

   print("volume averaging: start");
   timer_averaging->StartTimer();

   for (int i = 0; i < numberOfPoints; i++)
   {
      if (i%10000 == 0)
      {
         timer_inloop->StartTimer();
         print("point id = "+toString(i));
      }

      dataSet->GetPoint(i, center);
      pointLocator->FindPointsWithinRadius (radius, center, result);
      double numPoints = (double)result->GetNumberOfIds(); 

      for (int j = 0; j < (int)numPoints; j++)
      {
         vtkIdType id = result->GetId(j);

         Mvx[i]+=vxArray->GetTuple1(id);
         Mvy[i]+=vyArray->GetTuple1(id);
         Mvz[i]+=vzArray->GetTuple1(id);

         //std::cout << "j = "<<j<<std::endl;
      }

      if(numPoints > 0)
      {
         //mean values=integral value / domain size
         Mvx[i]=Mvx[i]/numPoints;
         Mvy[i]=Mvy[i]/numPoints;
         Mvz[i]=Mvz[i]/numPoints;
      }

      if (i%10000 == 0)
      {
         timer_inloop->StopTimer();
         print("numPoints = "+toString(numPoints));
         print("time: "+toString(timer_inloop->GetElapsedTime())+" s");
         print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
      }
   }

   timer_averaging->StopTimer();
   print("volume averaging: end");
   print("volume averaging time: "+toString(timer_averaging->GetElapsedTime())+" s");

   std::string vtkfilename = output;

   print("write data set to "+vtkfilename+": start");
   timer_write->StartTimer();

   vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridOut = vtkSmartPointer<vtkUnstructuredGrid>::New();
   unstructuredGridOut->CopyStructure(dataSet);

   vtkSmartPointer<vtkDoubleArray> mvxArray = vtkSmartPointer<vtkDoubleArray>::New();
   mvxArray->SetNumberOfComponents(1);
   mvxArray->SetName( "Mvx" );

   vtkSmartPointer<vtkDoubleArray> mvyArray = vtkSmartPointer<vtkDoubleArray>::New();
   mvyArray->SetNumberOfComponents(1);
   mvyArray->SetName( "Mvy" );

   vtkSmartPointer<vtkDoubleArray> mvzArray = vtkSmartPointer<vtkDoubleArray>::New();
   mvzArray->SetNumberOfComponents(1);
   mvzArray->SetName( "Mvz" );

   for(int i = 0; i < numberOfPoints; i++)
   {
      mvxArray->InsertNextValue(Mvx[i]);
      mvyArray->InsertNextValue(Mvy[i]);
      mvzArray->InsertNextValue(Mvz[i]);
   }

   unstructuredGridOut->GetPointData()->AddArray(mvxArray);
   unstructuredGridOut->GetPointData()->AddArray(mvyArray);
   unstructuredGridOut->GetPointData()->AddArray(mvzArray);

   vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
   writer->SetInputData(unstructuredGridOut);
   writer->SetFileName(vtkfilename.c_str());
   writer->SetDataModeToAscii();
   //writer->SetDataModeToBinary();
   writer->Update();

   print("write data set: end");
   timer_write->StopTimer();
   print("write data set time: "+toString(timer_write->GetElapsedTime())+" s");

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}


