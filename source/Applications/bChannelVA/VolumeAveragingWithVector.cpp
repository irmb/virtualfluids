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

//#define INDEX(x1, x2, x3, nx2, nx3) ( (nx2) * ( (nx3) * (x3) + (x2) ) + (x1) )

using namespace std;


vector<double> geo;

vector< vector <double> > zVx;
vector< vector <double> > zVy;
vector< vector <double> > zVz;
vector< vector <double> > zPress;

//surface averaged
vector<double> savVx;
vector<double> savVy;
vector<double> savVz;
vector<double> savPress;

//volume averaged
vector<double> vavVx;
vector<double> vavVy;
vector<double> vavVz;
vector<double> vavPress;

//fluctuations
//surface averaged
vector<double> savVxx;
vector<double> savVyy;
vector<double> savVzz;
vector<double> savPP;

//volume averaged
vector<double> vavVxx;
vector<double> vavVyy;
vector<double> vavVzz;
vector<double> vavPP;

//cross-correlations
vector<double> vavVxy;
vector<double> vavVxz;
vector<double> vavVyz;

//triple correlations
vector<double> vavVxyz;
vector<double> vavVxxx;
vector<double> vavVxxy;
vector<double> vavVxxz;
vector<double> vavVxyy;
vector<double> vavVxzz;
vector<double> vavVyyy;
vector<double> vavVyyz;
vector<double> vavVyzz;
vector<double> vavVzzz;

double deltax;
double origin[3];

double topo[] = {0.0, 10.0};

//////////////////////////////////////////////////////////////////////////
void LoadGeo(std::string dataName)
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

   for (double z = bounds[4]+deltax; z <= bounds[5]-deltax; z += deltax)
   {
      timer_calculate->StartTimer();
      print("calculate number of points: start");

      plane->SetOrigin(origin[0], origin[1], z);
      planeCut->Update();

      double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
      geo.push_back(numPoints);

      double org[3];
      plane->GetOrigin(org);
      print("origin = "+toString(org[0])+","+toString(org[1])+","+toString(org[2]));
      print("Z = "+toString(z));
      print("numPoints = "+toString(numPoints));
      print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
      timer_calculate->StopTimer();
      print("calculate number of points: end");
      print("calculate number of points time: "+toString(timer_calculate->GetElapsedTime())+" s");
   }

   int size = geo.size();

   zVx.resize(size);
   zVy.resize(size);
   zVz.resize(size);
   zPress.resize(size);

   savVx.resize(size);
   savVy.resize(size);
   savVz.resize(size);
   savPress.resize(size);

   vavVx.resize(size);
   vavVy.resize(size);
   vavVz.resize(size);
   vavPress.resize(size);

   savVxx.resize(size);
   savVyy.resize(size);
   savVzz.resize(size);
   savPP.resize(size);

   vavVxx.resize(size);
   vavVyy.resize(size);
   vavVzz.resize(size);
   vavPP.resize(size);

   vavVxy.resize(size);
   vavVxz.resize(size);
   vavVyz.resize(size);

   vavVxyz.resize(size);
   vavVxxx.resize(size);
   vavVxxy.resize(size);
   vavVxxz.resize(size);
   vavVxyy.resize(size);
   vavVxzz.resize(size);
   vavVyyy.resize(size);
   vavVyyz.resize(size);
   vavVyzz.resize(size);
   vavVzzz.resize(size);

   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void CalculateMeanMQVector(std::string dataName)
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

   int j = 0;

   for (double z = bounds[4]; z <= bounds[5]; z += deltax)
   {
      timer_calculate->StartTimer();
      print("calculate mean velocity: start");

      plane->SetOrigin(origin[0], origin[1], z);
      planeCut->Update();

      double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
      vtkSmartPointer<vtkDataArray> vxArray = planeCut->GetOutput()->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> pressArray = planeCut->GetOutput()->GetPointData()->GetArray("Press");
      double vxAvSumm = 0;
      double vyAvSumm = 0;
      double vzAvSumm = 0;
      double pressAvSumm = 0;

      for (int i = 0; i < numPoints; i++)
      {
         double vx = vxArray->GetTuple1(i);
         double vy = vyArray->GetTuple1(i);
         double vz = vzArray->GetTuple1(i);
         double press = pressArray->GetTuple1(i);

         vxAvSumm += vx;
         vyAvSumm += vy;
         vzAvSumm += vz;
         pressAvSumm += press;

         zVx[j].push_back(vx);
         zVy[j].push_back(vy);
         zVz[j].push_back(vz);
         zPress[j].push_back(press);
      }

      savVx[j] = vxAvSumm/geo[j];
      savVy[j] = vyAvSumm/geo[j];
      savVz[j] = vzAvSumm/geo[j];
      savPress[j] = pressAvSumm/geo[j];

      double org[3];
      plane->GetOrigin(org);
      print("origin = "+toString(org[0])+","+toString(org[1])+","+toString(org[2]));
      print("Z = "+toString(z));
      print("numPoints = "+toString(numPoints));
      print("vxAvSumm = "+toString(vxAvSumm));
      print("vxAv = "+toString(savVx[j]));
      print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
      timer_calculate->StopTimer();
      print("calculate mean velocity: end");
      print("calculate mean velocity time: "+toString(timer_calculate->GetElapsedTime())+" s");
      j++;
   }


   timer_global->StopTimer();
   print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
double G(double x, double l)
{
   if(fabs(x)<=l) 
      return l-fabs(x);
   else
      return 0.0;
}
//////////////////////////////////////////////////////////////////////////
double m(double x3, double l)
{
   return((G(x3,l))/(l*l));
}
//////////////////////////////////////////////////////////////////////////
void CorrectIndex(int& x, int nx)
{
   if (x < 0)   x = nx + x;
   if (x >= nx) x = x - nx;
}
//////////////////////////////////////////////////////////////////////////
void VolumeAveragingMQ(double L)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   print("volume averaging mq: start");
   timer_averaging->StartTimer();

   print("L = "+toString(L));

   int l = cint(L/deltax);

   print("l = "+toString(l));

   double mdv;

   int numberOfPoints = (int)savVx.size();

   for(int i = 0; i < numberOfPoints; i++)
   {
      double mvx   = 0;
      double mvy   = 0;
      double mvz   = 0;
      double mpress = 0;

      for(int z=-l;z<=+l;z++)
      {
         int zz = i+z;
         CorrectIndex(zz, numberOfPoints);

         mdv     = m((double)z, (double)l);
         mvx    += mdv*savVx[zz];
         mvy    += mdv*savVy[zz];
         mvz    += mdv*savVz[zz];
         mpress += mdv*savPress[zz];
      }

      vavVx[i] = mvx;
      vavVy[i] = mvy;
      vavVz[i] = mvz;
      vavPress[i] = mpress;
   }

   timer_averaging->StopTimer();
   print("volume averaging mq: end");
   print("volume averaging mq time: "+toString(timer_averaging->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void CalculateFluctuations()
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   print("calculate fluctuations: start");
   timer_averaging->StartTimer();

   int numberOfPoints = (int)vavVx.size();

   for (int i = 0; i < numberOfPoints; i++)
   {
      double vxxAvSumm = 0;
      double vyyAvSumm = 0;
      double vzzAvSumm = 0;
      double ppAvSumm = 0;

      int size = zVx[i].size();

      for (int j = 0; j < size; j++)
      {
         vxxAvSumm += zVx[i][j]-vavVx[i];
         vyyAvSumm += zVy[i][j]-vavVy[i];
         vzzAvSumm += zVz[i][j]-vavVz[i];
         ppAvSumm  += zPress[i][j]-vavPress[i];;
      }

      savVx[i] = vxxAvSumm/geo[i];
      savVy[i] = vyyAvSumm/geo[i];
      savVz[i] = vzzAvSumm/geo[i];
      savPP[i] =  ppAvSumm/geo[i];
   }

   timer_averaging->StopTimer();
   print("calculate fluctuations: end");
   print("calculate fluctuations time: "+toString(timer_averaging->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void VolumeAveragingF(double L)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
   print("volume averaging of fluctuations: start");
   timer_averaging->StartTimer();

   print("L = "+toString(L));

   int l = cint(L/deltax);

   print("l = "+toString(l));

   double mdv;

   int numberOfPoints = (int)savVx.size();

   for(int i = 0; i < numberOfPoints; i++)
   {
      double mvx   = 0;
      double mvy   = 0;
      double mvz   = 0;
      double mpress = 0;

      for(int z=-l;z<=+l;z++)
      {
         int zz = i+z;
         CorrectIndex(zz, numberOfPoints);

         mdv     = m((double)z, (double)l);
         mvx    += mdv*savVxx[zz];
         mvy    += mdv*savVyy[zz];
         mvz    += mdv*savVzz[zz];
         mpress += mdv*savPP[zz];
      }

      vavVxx[i] = mvx;
      vavVyy[i] = mvy;
      vavVzz[i] = mvz;
      vavPP[i]  = mpress;
   }

   timer_averaging->StopTimer();
   print("volume averaging of fluctuations: end");
   print("volume averaging of fluctuations time: "+toString(timer_averaging->GetElapsedTime())+" s");
}
//////////////////////////////////////////////////////////////////////////
void CalculateCorrelations()
{
   int numberOfPoints = (int)vavVxx.size();

   for(int i = 0; i < numberOfPoints; i++)
   {
      vavVxy[i] = vavVxx[i] * vavVyy[i];
      vavVxz[i] = vavVxx[i] * vavVzz[i];
      vavVyz[i] = vavVyy[i] * vavVzz[i];

      vavVxyz[i] = vavVxx[i] * vavVyy[i] * vavVzz[i];
      vavVxxx[i] = vavVxx[i] * vavVxx[i] * vavVxx[i];
      vavVxxy[i] = vavVxx[i] * vavVxx[i] * vavVyy[i];
      vavVxxz[i] = vavVxx[i] * vavVxx[i] * vavVzz[i];
      vavVxyy[i] = vavVxx[i] * vavVyy[i] * vavVyy[i];
      vavVxzz[i] = vavVxx[i] * vavVzz[i] * vavVzz[i];
      vavVyyy[i] = vavVyy[i] * vavVyy[i] * vavVyy[i];
      vavVyyz[i] = vavVyy[i] * vavVyy[i] * vavVzz[i];
      vavVyzz[i] = vavVxx[i] * vavVzz[i] * vavVzz[i];
      vavVzzz[i] = vavVzz[i] * vavVzz[i] * vavVzz[i];
   }
}
//////////////////////////////////////////////////////////////////////////
void WriteData(std::string output)
{
   FILE * pFile;
   std::string fname = output+".dat";

   pFile = fopen (fname.c_str(),"w");
   if(pFile == NULL)
   {
      printf("Can not open file %s!", fname.c_str());
      return;
   }

   fprintf (pFile, "z;avVx;avVy;avVz;avPress;avVxx;avVyy;avVzz;avPP;avVxy;avVxz;avVyz;avVxyz;avVxxx;avVxxy;avVxxz;avVxyy;avVxzz;avVyyy;avVyyz;avVyzz;avVzzz\n");

   int size = vavVx.size();

   for (int i = 0; i < size; i++)
   {
      double z =origin[2] + (double)i*deltax;
      fprintf(pFile, "%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n",z,geo[i],vavVx[i],vavVy[i],vavVz[i],vavPress[i],vavVxx[i],vavVyy[i],vavVzz[i],vavPP[i],vavVxy[i],vavVxz[i],vavVyz[i],vavVxyz[i],vavVxxx[i],vavVxxy[i],vavVxxz[i],vavVxyy[i],vavVxzz[i],vavVyyy[i],vavVyyz[i],vavVyzz[i],vavVzzz[i]); 
   }


   fclose (pFile);                                                                    
}                                                                                     
//////////////////////////////////////////////////////////////////////////
void VolumeAveragingWithVector()
{
   deltax = 2.66223; //cdx
   //deltax = 0.665557; //fdx

   LoadGeo("/work/koskuche/SFB880/BKanalTestBerlin/grid/nodes_0.pvtu");

   CalculateMeanMQVector("/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/steps/step_488750.pvtu");

   CalculateFluctuations();

   CalculateCorrelations();

   WriteData("/work/koskuche/SFB880/BKanalTestBerlin/avData/avData1.csv");
}

