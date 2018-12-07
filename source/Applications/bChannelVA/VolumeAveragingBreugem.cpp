#include "Postprocessing.h"
#include <vector>
#include <vtkDataArray.h>
#include <vtkPointLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkTimerLog.h>
#include "MemoryUtil.h"

//////////////////////////////////////////////////////////////////////////
//df = 10 m, dp = 10.02 m
//////////////////////////////////////////////////////////////////////////
void VolumeAveragingBreugem(std::string dataName, double L, double dx, std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_locator = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();
   timer_global->StartTimer();

   //double l=2.0*(df+dp);

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

   double bounds[6];
   dataSet->GetBounds(bounds);

   double X_Min = bounds[0];
   double X_Max = bounds[1];
   double Y_Min = bounds[2];
   double Y_Max = bounds[3];
   double Z_Min = bounds[4];
   double Z_Max = bounds[5];

   Mvx.resize(numberOfPoints,0.0);
   Mvy.resize(numberOfPoints,0.0);
   Mvz.resize(numberOfPoints,0.0);


   vtkSmartPointer<vtkDataArray> vxArray = dataSet->GetPointData()->GetArray("Mvx");
   vtkSmartPointer<vtkDataArray> vyArray = dataSet->GetPointData()->GetArray("Mvy");
   vtkSmartPointer<vtkDataArray> vzArray = dataSet->GetPointData()->GetArray("Mvz");

   double X[3];
   
   double x_min, x_max,y_min,y_max,z_min,z_max;
   double dy,dz;
   double l = L/2.0;

   x_min=x_max=L/2.0;
   z_min=z_max=L/2.0;
   y_min=L/2.0;
   y_max=L/2.0;
   //double lx_mean=(x_min+x_max)/2.0;
   //double ly_mean=(y_min+y_max)/2.0;
   //double lz_mean=(z_min+z_max)/2.0;
   dy=dz=dx; 

   double dx3 = dx*dx*dx;
   //double norm = (lx_mean*lx_mean)*(ly_mean*ly_mean)*(lz_mean*lz_mean);

   print("volume averaging: start");
   timer_averaging->StartTimer();

   for (int i = 0; i < numberOfPoints; i++)
   {
      if (i%10000 == 0)
      {
         timer_inloop->StartTimer();
         print("point id = "+toString(i));
      }

      dataSet->GetPoint(i, X);

      if(X[0] < X_Min+l || X[0] > X_Max-l || X[1] < Y_Min+l || X[1] > Y_Max-l || X[2] < Z_Min+l || X[2] > Z_Max-l)
         continue;

      //double temp;
      //if (temp=fabs(X[0]-X_Min) < l)
      //{
      //   x_min = temp;
      //}
      //else
      //   x_min = l;
      //if (temp=fabs(X[1]-Y_Min) < l)
      //{
      //   y_min = temp;
      //}
      //else
      //   y_min = l;
      //if (temp=fabs(X[2]-Z_Min) < l)
      //{
      //   z_min = temp;
      //}
      //else
      //   z_min = l;
     
      //if (temp=fabs(X[0]-X_Max) < l)
      //{
      //   x_max = temp;
      //}
      //else
      //   x_max = l;
      //if (temp=fabs(X[1]-Y_Max) < l)
      //{
      //   y_max = temp;
      //}
      //else
      //   y_max = l;
      //if (temp=fabs(X[2]-Z_Max) < l)
      //{
      //   z_max = temp;
      //}
      //else
      //   z_max = l;

      double lx_mean=(x_min+x_max)/2.0;
      double ly_mean=(y_min+y_max)/2.0;
      double lz_mean=(z_min+z_max)/2.0;

      double norm = (lx_mean*lx_mean)*(ly_mean*ly_mean)*(lz_mean*lz_mean);

      int numPoints = 0;

      for(double x=X[0]-x_min;x<X[0]+x_max;x+=dx){
         for(double y=X[1]-y_min;y<X[1]+y_max;y+=dy){
            for(double z=X[2]-z_min;z<X[2]+z_max;z+=dz){
               vtkIdType id = dataSet->FindPoint(X[0]+x, X[1]+y, X[2]+z);
               if( id >= 0)
               {
                  double mdv =m(x+dx/2.0-X[0],y+dy/2.0-X[1],z+dz/2.0-X[2],x_min,x_max,y_min,y_max,z_min,z_max,norm)*dx3;
                  Mvx[i]+=mdv*vxArray->GetTuple1(id);
                  Mvy[i]+=mdv*vyArray->GetTuple1(id);
                  Mvz[i]+=mdv*vzArray->GetTuple1(id);
                  numPoints++;
               }

            }
         }
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

   std::string vtkfilename = output+".vti";

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


