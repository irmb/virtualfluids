#include "Postprocessing.h"
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkTimerLog.h>
#include "MemoryUtil.h"

void SurfaceAveragingMeanVelocity(std::string dataName, double step, std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
   timer_global->StartTimer();

   timer_read->StartTimer();
   print("read data set from "+dataName+": start");
   vtkSmartPointer<vtkDataSet> dataSet(ReadDataSet(dataName.c_str()));
   print("read data set from "+dataName+": end");
   timer_read->StopTimer();
   print("read data set time: "+toString(timer_read->GetElapsedTime())+" s");

   double bounds[6];
   dataSet->GetBounds(bounds);
   double center[3];
   dataSet->GetCenter(center);

   double origin[3];
   origin[0] = center[0];
   origin[1] = center[1];
   origin[2] = center[2];

   //vtkPlane
   vtkSmartPointer<vtkPlane> plane = vtkPlane::New();
   plane->SetNormal(0.0, 0.0, 1.0);

   //Cut
   vtkSmartPointer<vtkCutter> planeCut = vtkCutter::New();
   planeCut->SetInputData(dataSet);
   planeCut->SetCutFunction(plane);

   std::ofstream ostr;
   std::string fname = output+".dat";
   ostr.open(fname.c_str(), std::ios_base::out);
   if(!ostr)
   { 
      ostr.clear();
   }
   ostr << "z" << "\t";
   ostr << "<Vx>" << std::endl;

   for (double j = bounds[4]; j <= bounds[5]; j += step)
   {
      timer_calculate->StartTimer();
      print("calculate mean velocity: start");

      plane->SetOrigin(origin[0], origin[1], j);
      planeCut->Update();

      double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
      vtkSmartPointer<vtkDataArray> vxArray = planeCut->GetOutput()->GetPointData()->GetArray("Mvx");
      double vxAvSumm = 0;
      double vxAv = 0;
      for (int i = 0; i < numPoints; i++)
      {
         vxAvSumm += vxArray->GetTuple1(i);
      }

      if(numPoints > 0)
      {
         vxAv = vxAvSumm/numPoints;
      }

      ostr << j << "\t";
      ostr << vxAv << std::endl;

      double org[3];
      plane->GetOrigin(org);
      print("origin = "+toString(org[0])+","+toString(org[1])+","+toString(org[2]));
      print("Z = "+toString(j));
      print("numPoints = "+toString(numPoints));
      print("vxAvSumm = "+toString(vxAvSumm));
      print("vxAv = "+toString(vxAv));
      print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
      timer_calculate->StopTimer();
      print("calculate mean velocity: end");
      print("calculate mean velocity time: "+toString(timer_calculate->GetElapsedTime())+" s");
   }

   ostr.close();

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}

