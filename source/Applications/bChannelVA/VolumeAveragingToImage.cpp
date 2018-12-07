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
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include "MemoryUtil.h"

void VolumeAveragingToImage(std::string dataName, double L, double dx, double image_dx, std::string output)
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
   
   print("read data set from "+dataName+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");
   
   double l = L/2.0;

   vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
   double delta = image_dx;
   double bounds[6];
   dataSet->GetBounds(bounds);
   int extent[6];
   extent[0] = 0;
   extent[2] = 0;
   extent[4] = 0;
   extent[1] = (int)(fabs((bounds[1] - bounds[0]-L)/delta));
   extent[3] = (int)(fabs((bounds[3] - bounds[2]-L)/delta));
   extent[5] = (int)(fabs((bounds[5] - bounds[4]-L)/delta));
   image->SetExtent(extent);
   image->SetOrigin(bounds[0]+l, bounds[2]+l, bounds[4]+l);
   image->SetSpacing(delta, delta, delta);

   std::string vtkfilename = output+"_s.vti";
   vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
   writer->SetInputData(image);
   writer->SetFileName(vtkfilename.c_str());
   writer->SetDataModeToAscii();
   writer->Update();

   int numberOfPoints = image->GetNumberOfPoints();

   print("write image structure to "+vtkfilename);
   print("number of points in image = "+toString(numberOfPoints));

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

   double X[3];

   double x_min, x_max,y_min,y_max,z_min,z_max;
   double dy,dz;
   

   x_min=x_max=L/2.0;
   z_min=z_max=L/2.0;
   y_min=L/2.0;
   y_max=L/2.0;
   double lx_mean=(x_min+x_max)/2.0;
   double ly_mean=(y_min+y_max)/2.0;
   double lz_mean=(z_min+z_max)/2.0;
   dy=dz=dx; 
   double dx3 = dx*dx*dx;
   double norm = (lx_mean*lx_mean)*(ly_mean*ly_mean)*(lz_mean*lz_mean);

   print("volume averaging: start");
   timer_averaging->StartTimer();

   int numPoints = 0;
   double coord[3];

   for (int i = 0; i < numberOfPoints; i++)
   {
      if (i%100 == 0)
      {
         print("point id = "+toString(i));
         timer_inloop->StartTimer();
         
      }

      image->GetPoint(i, X);

      for(double z=X[2]-z_min;z<X[2]+z_max;z+=dz){
         for(double y=X[1]-y_min;y<X[1]+y_max;y+=dy){
            for(double x=X[0]-x_min;x<X[0]+x_max;x+=dx){   
               //vtkIdType id = dataSet->FindPoint(X[0]+x, X[1]+y, X[2]+z);
               coord[0]=X[0]+x;
               coord[1]=X[1]+y;
               coord[2]=X[2]+z;
               vtkIdType id = pointLocator->FindClosestPoint(coord);
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

      if (i%100 == 0)
      {
         timer_inloop->StopTimer();
         print("numPoints = "+toString(numPoints));
         print("time: "+toString(timer_inloop->GetElapsedTime())+" s");
         print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe()/1e9)+" GByte");
         numPoints = 0;
      }
   }

   timer_averaging->StopTimer();
   print("volume averaging: end");
   print("volume averaging time: "+toString(timer_averaging->GetElapsedTime())+" s");

   vtkfilename = output+"_f.vti";

   print("write data set to "+vtkfilename+": start");
   timer_write->StartTimer();

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

   image->GetPointData()->AddArray(mvxArray);
   image->GetPointData()->AddArray(mvyArray);
   image->GetPointData()->AddArray(mvzArray);

   //vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
   //writer->SetInputData(image);
   writer->SetFileName(vtkfilename.c_str());
   //writer->SetDataModeToAscii();
   writer->Update();

   print("write data set: end");
   timer_write->StopTimer();
   print("write data set time: "+toString(timer_write->GetElapsedTime())+" s");

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}


