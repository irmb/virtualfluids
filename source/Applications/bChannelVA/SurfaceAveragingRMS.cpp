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

void SurfaceAveragingRMS(std::string dataName, double step, std::string output)
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
   ostr << "<Vxx>" << "\t";
   ostr << "<Vyy>" << "\t";
   ostr << "<Vzz>" << std::endl;

   for (double j = bounds[4]; j <= bounds[5]; j += step)
   {
      timer_calculate->StartTimer();
      print("calculate mean velocity: start");

      plane->SetOrigin(origin[0], origin[1], j);
      planeCut->Update();

      double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
      vtkSmartPointer<vtkDataArray> vxxArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxx");
      double vxxAvSumm = 0;
      double vxxAv = 0;
      vtkSmartPointer<vtkDataArray> vyyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vyy");
      double vyyAvSumm = 0;
      double vyyAv = 0;
      vtkSmartPointer<vtkDataArray> vzzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vzz");
      double vzzAvSumm = 0;
      double vzzAv = 0;
      for (int i = 0; i < numPoints; i++)
      {
         vxxAvSumm += vxxArray->GetTuple1(i);
         vyyAvSumm += vyyArray->GetTuple1(i);
         vzzAvSumm += vzzArray->GetTuple1(i);
      }

      if(numPoints > 0)
      {
         vxxAv = vxxAvSumm/numPoints;
         vyyAv = vyyAvSumm/numPoints;
         vzzAv = vzzAvSumm/numPoints;
      }

      ostr << j << "\t";
      ostr << vxxAv << "\t";
      ostr << vyyAv << "\t";
      ostr << vzzAv << std::endl;

      double org[3];
      plane->GetOrigin(org);
      print("origin = "+toString(org[0])+","+toString(org[1])+","+toString(org[2]));
      print("Z = "+toString(j));
      print("numPoints = "+toString(numPoints));
      print("vxxAvSumm = "+toString(vxxAvSumm));
      print("vxxAv = "+toString(vxxAv));
      print("vyyAvSumm = "+toString(vyyAvSumm));
      print("vyyAv = "+toString(vyyAv));
      print("vzzAvSumm = "+toString(vzzAvSumm));
      print("vzzAv = "+toString(vzzAv));
      print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe()/1e9)+" GByte");
      timer_calculate->StopTimer();
      print("calculate mean velocity: end");
      print("calculate mean velocity time: "+toString(timer_calculate->GetElapsedTime())+" s");
   }

   ostr.close();

   print("write data to "+fname);
   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}

