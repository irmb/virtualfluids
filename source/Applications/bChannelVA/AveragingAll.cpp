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

using namespace std;

namespace AA
{

   //vector<double> z;
   //vector<double> vx;
   //vector<double> vy;
   //vector<double> vz;

   void SurfaceAveragingAll(std::string dataName, double start, double step, std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
      timer_global->StartTimer();

      timer_read->StartTimer();
      print("read data set from "+dataName+": start");
      vtkDataSet *dataSet = ReadDataSet(dataName.c_str());
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
      vtkPlane *plane = vtkPlane::New();
      plane->SetNormal(0.0, 0.0, 1.0);

      //Cut
      vtkCutter *planeCut = vtkCutter::New();
      planeCut->SetInputData(dataSet);
      planeCut->SetCutFunction(plane);

      std::string fname = output;

      FILE *file;
      file = fopen(fname.c_str(), "w");

      if (file == NULL)
      {
         print("can not open " + fname);
         return;
      }

      fprintf(file, "z;<Vx>;<Vy>;<Vz>;<Vxx>;<Vyy>;<Vzz>;<Vxy>;<Vxz>;<Vyz>\n");

      //for (double j = bounds[4]; j <= bounds[5]; j += step)
      for (double j = start; j <= bounds[5]; j += step)
      {
         //timer_calculate->StartTimer();
         //print("calculate plane: start");

         plane->SetOrigin(origin[0], origin[1], j);
         planeCut->Update();

         double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
         if (numPoints < 1)
         {
            continue;
         }
         vtkDataArray *vxArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVx");
         vtkDataArray *vyArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVy");
         vtkDataArray *vzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVz");
         vtkDataArray *vxxArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVxx");
         vtkDataArray *vyyArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVyy");
         vtkDataArray *vzzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVzz");
         vtkDataArray *vxyArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVxy");
         vtkDataArray *vxzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVxz");
         vtkDataArray *vyzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVyz");
         double vxAvSumm = 0;
         double vxAv = 0;
         double vyAvSumm = 0;
         double vyAv = 0;
         double vzAvSumm = 0;
         double vzAv = 0;
         double vxxAvSumm = 0;
         double vxxAv = 0;
         double vyyAvSumm = 0;
         double vyyAv = 0;
         double vzzAvSumm = 0;
         double vzzAv = 0;
         double vxyAvSumm = 0;
         double vxyAv = 0;
         double vxzAvSumm = 0;
         double vxzAv = 0;
         double vyzAvSumm = 0;
         double vyzAv = 0;

         for (int i = 0; i < numPoints; i++)
         {
            vxAvSumm += vxArray->GetTuple1(i);
            vyAvSumm += vyArray->GetTuple1(i);
            vzAvSumm += vzArray->GetTuple1(i);
            vxxAvSumm += vxxArray->GetTuple1(i);
            vyyAvSumm += vyyArray->GetTuple1(i);
            vzzAvSumm += vzzArray->GetTuple1(i);
            vxyAvSumm += vxyArray->GetTuple1(i);
            vxzAvSumm += vxzArray->GetTuple1(i);
            vyzAvSumm += vyzArray->GetTuple1(i);
         }

         if (numPoints > 0)
         {
            vxAv = vxAvSumm/numPoints;
            vyAv = vyAvSumm/numPoints;
            vzAv = vzAvSumm/numPoints;
            vxxAv = vxxAvSumm/numPoints;
            vyyAv = vyyAvSumm/numPoints;
            vzzAv = vzzAvSumm/numPoints;
            vxyAv = vxyAvSumm/numPoints;
            vxzAv = vxzAvSumm/numPoints;
            vyzAv = vyzAvSumm/numPoints;
         }

         fprintf(file, "%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n", j, vxAv, vyAv, vzAv, vxxAv, vyyAv, vzzAv, vxyAv, vxzAv, vyzAv);

         //double org[3];
         //plane->GetOrigin(org);
         //print("origin = "+toString(org[0])+","+toString(org[1])+","+toString(org[2]));
         //print("Z = "+toString(j));
         //print("numPoints = "+toString(numPoints));
         //print("vxAvSumm = "+toString(vxAvSumm));
         //print("vxAv = "+toString(vxAv));
         //print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
         //timer_calculate->StopTimer();
         //print("calculate plane: end");
         //print("calculate plane time: "+toString(timer_calculate->GetElapsedTime())+" s");
      }

      //ostr.close();
      fclose(file);

      planeCut->Delete();
      plane->Delete();
      dataSet->Delete();

      timer_global->StopTimer();
      print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
   }
   //////////////////////////////////////////////////////////////////////////
   void SurfaceAveragingAllWithLevels(std::string dataName, int minLevel, int maxLevel, vector<double> levelCoords, double corseDeltax, std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
      timer_global->StartTimer();

      timer_read->StartTimer();
      print("read data set from "+dataName+": start");
      vtkDataSet* dataSet = ReadDataSet(dataName.c_str());
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
      vtkPlane* plane = vtkPlane::New();
      plane->SetNormal(0.0, 0.0, 1.0);

      //Cut
      vtkCutter* planeCut = vtkCutter::New();
      planeCut->SetInputData(dataSet);
      planeCut->SetCutFunction(plane);

      std::string fname = output;

      FILE *file;
      file = fopen(fname.c_str(), "w");

      if (file == NULL)
      {
         print("can not open " + fname);
         return;
      }

      fprintf(file, "z;<Vx>;<Vy>;<Vz>;<Vxx>;<Vyy>;<Vzz>;<Vxy>;<Vxz>;<Vyz>\n");

      for (int level = maxLevel; level >= minLevel; level--)
      {
         double step = corseDeltax/(double)(1<<level);
         double start = levelCoords[level]+step*0.5;
         double stop;
         if (level > minLevel)
         {
            stop = levelCoords[level-1] - step*0.5;
         }
         else
         {
            stop = bounds[5];
         }

         timer_calculate->StartTimer();
         print("calculate level "+toString(level)+" : start");
         print("step = "+toString(step));


         for (double j = start; j <=stop; j += step)
         {
            //timer_calculate->StartTimer();
            //print("calculate plane: start");

            plane->SetOrigin(origin[0], origin[1], j);
            planeCut->Update();

            double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
            if (numPoints < 1)
            {
               continue;
            }
            vtkSmartPointer<vtkDataArray> vxArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVx");
            vtkSmartPointer<vtkDataArray> vyArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVy");
            vtkSmartPointer<vtkDataArray> vzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVz");
            vtkSmartPointer<vtkDataArray> vxxArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVxx");
            vtkSmartPointer<vtkDataArray> vyyArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVyy");
            vtkSmartPointer<vtkDataArray> vzzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVzz");
            vtkSmartPointer<vtkDataArray> vxyArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVxy");
            vtkSmartPointer<vtkDataArray> vxzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVxz");
            vtkSmartPointer<vtkDataArray> vyzArray = planeCut->GetOutput()->GetPointData()->GetArray("AvVyz");
            double vxAvSumm = 0;
            double vxAv = 0;
            double vyAvSumm = 0;
            double vyAv = 0;
            double vzAvSumm = 0;
            double vzAv = 0;
            double vxxAvSumm = 0;
            double vxxAv = 0;
            double vyyAvSumm = 0;
            double vyyAv = 0;
            double vzzAvSumm = 0;
            double vzzAv = 0;
            double vxyAvSumm = 0;
            double vxyAv = 0;
            double vxzAvSumm = 0;
            double vxzAv = 0;
            double vyzAvSumm = 0;
            double vyzAv = 0;

            for (int i = 0; i < numPoints; i++)
            {
               vxAvSumm += vxArray->GetTuple1(i);
               vyAvSumm += vyArray->GetTuple1(i);
               vzAvSumm += vzArray->GetTuple1(i);
               vxxAvSumm += vxxArray->GetTuple1(i);
               vyyAvSumm += vyyArray->GetTuple1(i);
               vzzAvSumm += vzzArray->GetTuple1(i);
               vxyAvSumm += vxyArray->GetTuple1(i);
               vxzAvSumm += vxzArray->GetTuple1(i);
               vyzAvSumm += vyzArray->GetTuple1(i);
            }

            if (numPoints > 0)
            {
               vxAv = vxAvSumm/numPoints;
               vyAv = vyAvSumm/numPoints;
               vzAv = vzAvSumm/numPoints;
               vxxAv = vxxAvSumm/numPoints;
               vyyAv = vyyAvSumm/numPoints;
               vzzAv = vzzAvSumm/numPoints;
               vxyAv = vxyAvSumm/numPoints;
               vxzAv = vxzAvSumm/numPoints;
               vyzAv = vyzAvSumm/numPoints;
            }

            fprintf(file, "%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n", j, vxAv, vyAv, vzAv, vxxAv, vyyAv, vzzAv, vxyAv, vxzAv, vyzAv);

            //double org[3];
            //plane->GetOrigin(org);
            //print("origin = "+toString(org[0])+","+toString(org[1])+","+toString(org[2]));
            //print("Z = "+toString(j));
            //print("numPoints = "+toString(numPoints));
            //print("vxAvSumm = "+toString(vxAvSumm));
            //print("vxAv = "+toString(vxAv));
            //print("actual memory usage: "+toString(MemoryUtil::getPhysMemUsedByMe())+" Byte");
            //timer_calculate->StopTimer();
            //print("calculate plane: end");
            //print("calculate plane time: "+toString(timer_calculate->GetElapsedTime())+" s");
         }
         print("calculate level "+toString(level)+" : stop");
         timer_calculate->StopTimer();
         print("calculate level time: "+toString(timer_calculate->GetElapsedTime())+" s");
      }

      fclose(file);

      planeCut->Delete();
      plane->Delete();
      dataSet->Delete();

      timer_global->StopTimer();
      print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
   }
   //////////////////////////////////////////////////////////////////////////
   void SurfaceAveragingVelocityWithLevels(std::string dataName, int minLevel, int maxLevel, vector<double> levelCoords, double corseDeltax, std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
      timer_global->StartTimer();

      timer_read->StartTimer();
      print("read data set from "+dataName+": start");
      vtkDataSet* dataSet = ReadDataSet(dataName.c_str());
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
      vtkPlane* plane = vtkPlane::New();
      plane->SetNormal(0.0, 0.0, 1.0);

      //Cut
      vtkCutter* planeCut = vtkCutter::New();
      planeCut->SetInputData(dataSet);
      planeCut->SetCutFunction(plane);

      std::string fname = output;

      FILE *file;
      file = fopen(fname.c_str(), "w");

      if (file == NULL)
      {
         print("can not open " + fname);
         return;
      }

      fprintf(file, "z;<Vx>;<Vy>;<Vz>\n");

      for (int level = maxLevel; level >= minLevel; level--)
      {
         double step = corseDeltax/(double)(1<<level);
         double start = levelCoords[level]+step*0.5;
         double stop;
         if (level > minLevel)
         {
            stop = levelCoords[level-1] - step*0.5;
         }
         else
         {
            stop = bounds[5];
         }

         timer_calculate->StartTimer();
         print("calculate level "+toString(level)+" : start");
         print("step = "+toString(step));


         for (double j = start; j <=stop; j += step)
         {
            plane->SetOrigin(origin[0], origin[1], j);
            planeCut->Update();

            double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
            if (numPoints < 1)
            {
               continue;
            }
            vtkDataArray *vxArray = planeCut->GetOutput()->GetPointData()->GetArray("Mvx");
            vtkDataArray *vyArray = planeCut->GetOutput()->GetPointData()->GetArray("Mvy");
            vtkDataArray *vzArray = planeCut->GetOutput()->GetPointData()->GetArray("Mvz");

            double vxAvSumm = 0;
            double vxAv = 0;
            double vyAvSumm = 0;
            double vyAv = 0;
            double vzAvSumm = 0;
            double vzAv = 0;

            for (int i = 0; i < numPoints; i++)
            {
               vxAvSumm += vxArray->GetTuple1(i);
               vyAvSumm += vyArray->GetTuple1(i);
               vzAvSumm += vzArray->GetTuple1(i);
            }

            if (numPoints > 0)
            {
               vxAv = vxAvSumm/numPoints;
               vyAv = vyAvSumm/numPoints;
               vzAv = vzAvSumm/numPoints;
            }

            fprintf(file, "%g;%g;%g;%g\n", j, vxAv, vyAv, vzAv);

         }
         print("calculate level "+toString(level)+" : stop");
         timer_calculate->StopTimer();
         print("calculate level time: "+toString(timer_calculate->GetElapsedTime())+" s");
      }

      fclose(file);

      planeCut->Delete();
      plane->Delete();
      dataSet->Delete();

      timer_global->StopTimer();
      print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
   }
   //////////////////////////////////////////////////////////////////////////
   void SurfaceAveragingFluctuationsWithLevels(std::string dataName, int minLevel, int maxLevel, vector<double> levelCoords, double corseDeltax, std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_global = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_read = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_calculate = vtkSmartPointer<vtkTimerLog>::New();
      timer_global->StartTimer();

      timer_read->StartTimer();
      print("read data set from "+dataName+": start");
      vtkDataSet* dataSet = ReadDataSet(dataName.c_str());
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
      vtkPlane* plane = vtkPlane::New();
      plane->SetNormal(0.0, 0.0, 1.0);

      //Cut
      vtkCutter* planeCut = vtkCutter::New();
      planeCut->SetInputData(dataSet);
      planeCut->SetCutFunction(plane);

      std::string fname = output;

      FILE *file;
      file = fopen(fname.c_str(), "w");

      if (file == NULL)
      {
         print("can not open " + fname);
         return;
      }

      fprintf(file, "z;Vxx;Vyy;Vzz;Vxy;Vxz;Vyz;Vxxx;Vxxy;Vxxz;Vyyy;Vyyx;Vyyz;Vzzz;Vzzx;Vzzy;Vxyz\n");

      for (int level = maxLevel; level >= minLevel; level--)
      {
         double step = corseDeltax/(double)(1<<level);
         double start = levelCoords[level]+step*0.5;
         double stop;
         if (level > minLevel)
         {
            stop = levelCoords[level-1] - step*0.5;
         }
         else
         {
            stop = bounds[5];
         }

         timer_calculate->StartTimer();
         print("calculate level "+toString(level)+" : start");
         print("step = "+toString(step));


         for (double j = start; j <=stop; j += step)
         {
            //timer_calculate->StartTimer();
            //print("calculate plane: start");

            plane->SetOrigin(origin[0], origin[1], j);
            planeCut->Update();

            double numPoints = planeCut->GetOutput()->GetNumberOfPoints();
            if (numPoints < 1)
            {
               continue;
            }

            vtkSmartPointer<vtkDataArray> vxxArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxx");
            vtkSmartPointer<vtkDataArray> vyyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vyy");
            vtkSmartPointer<vtkDataArray> vzzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vzz");
            vtkSmartPointer<vtkDataArray> vxyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxy");
            vtkSmartPointer<vtkDataArray> vxzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxz");
            vtkSmartPointer<vtkDataArray> vyzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vyz");

            vtkSmartPointer<vtkDataArray> vxxxArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxxx");
            vtkSmartPointer<vtkDataArray> vxxyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxxy");
            vtkSmartPointer<vtkDataArray> vxxzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxxz");
            vtkSmartPointer<vtkDataArray> vyyyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vyyy");
            vtkSmartPointer<vtkDataArray> vyyxArray = planeCut->GetOutput()->GetPointData()->GetArray("Vyyx");
            vtkSmartPointer<vtkDataArray> vyyzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vyyz");
            vtkSmartPointer<vtkDataArray> vzzzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vzzz");
            vtkSmartPointer<vtkDataArray> vzzxArray = planeCut->GetOutput()->GetPointData()->GetArray("Vzzx");
            vtkSmartPointer<vtkDataArray> vzzyArray = planeCut->GetOutput()->GetPointData()->GetArray("Vzzy");
            vtkSmartPointer<vtkDataArray> vxyzArray = planeCut->GetOutput()->GetPointData()->GetArray("Vxyz");

            double vxxAvSumm = 0;
            double vyyAvSumm = 0;
            double vzzAvSumm = 0;
            double vxyAvSumm = 0;
            double vxzAvSumm = 0;
            double vyzAvSumm = 0;

            double vxxxAvSumm = 0;
            double vxxyAvSumm = 0;
            double vxxzAvSumm = 0;
            double vyyyAvSumm = 0;
            double vyyxAvSumm = 0;
            double vyyzAvSumm = 0;
            double vzzzAvSumm = 0;
            double vzzxAvSumm = 0;
            double vzzyAvSumm = 0;
            double vxyzAvSumm = 0;

            double vxxAv = 0;
            double vyyAv = 0;
            double vzzAv = 0;
            double vxyAv = 0;
            double vxzAv = 0;
            double vyzAv = 0;

            double vxxxAv = 0;
            double vxxyAv = 0;
            double vxxzAv = 0;
            double vyyyAv = 0;
            double vyyxAv = 0;
            double vyyzAv = 0;
            double vzzzAv = 0;
            double vzzxAv = 0;
            double vzzyAv = 0;
            double vxyzAv = 0;


            for (int i = 0; i < numPoints; i++)
            {
               vxxAvSumm += vxxArray->GetTuple1(i);
               vyyAvSumm += vyyArray->GetTuple1(i);
               vzzAvSumm += vzzArray->GetTuple1(i);
               vxyAvSumm += vxyArray->GetTuple1(i);
               vxzAvSumm += vxzArray->GetTuple1(i);
               vyzAvSumm += vyzArray->GetTuple1(i);

               vxxxAvSumm += vxxxArray->GetTuple1(i);
               vxxyAvSumm += vxxyArray->GetTuple1(i);
               vxxzAvSumm += vxxzArray->GetTuple1(i);
               vyyyAvSumm += vyyyArray->GetTuple1(i);
               vyyxAvSumm += vyyxArray->GetTuple1(i);
               vyyzAvSumm += vyyzArray->GetTuple1(i);
               vzzzAvSumm += vzzzArray->GetTuple1(i);
               vzzxAvSumm += vzzxArray->GetTuple1(i);
               vzzyAvSumm += vzzyArray->GetTuple1(i);
               vxyzAvSumm += vxyzArray->GetTuple1(i);
            }

            if (numPoints > 0)
            {
               vxxAv = vxxAvSumm/numPoints;
               vyyAv = vyyAvSumm/numPoints;
               vzzAv = vzzAvSumm/numPoints;
               vxyAv = vxyAvSumm/numPoints;
               vxzAv = vxzAvSumm/numPoints;
               vyzAv = vyzAvSumm/numPoints;

               vxxxAv = vxxxAvSumm/numPoints;
               vxxyAv = vxxyAvSumm/numPoints;
               vxxzAv = vxxzAvSumm/numPoints;
               vyyyAv = vyyyAvSumm/numPoints;
               vyyxAv = vyyxAvSumm/numPoints;
               vyyzAv = vyyzAvSumm/numPoints;
               vzzzAv = vzzzAvSumm/numPoints;
               vzzxAv = vzzxAvSumm/numPoints;
               vzzyAv = vzzyAvSumm/numPoints;
               vxyzAv = vxyzAvSumm/numPoints;
            }

            fprintf(file, "%g;", j);
            fprintf(file, "%g;%g;%g;%g;%g;%g;", vxxAv, vyyAv, vzzAv, vxyAv, vxzAv, vyzAv);
            fprintf(file, "%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n", vxxxAv, vxxyAv, vxxzAv, vyyyAv, vyyxAv, vyyzAv, vzzzAv, vzzxAv, vzzyAv, vxyzAv);

         }
         print("calculate level "+toString(level)+" : stop");
         timer_calculate->StopTimer();
         print("calculate level time: "+toString(timer_calculate->GetElapsedTime())+" s");
      }

      fclose(file);

      planeCut->Delete();
      plane->Delete();
      dataSet->Delete();

      timer_global->StopTimer();
      print("Total time: "+toString(timer_global->GetElapsedTime())+" s");
   }
   //////////////////////////////////////////////////////////////////////////
   //void ReadCSV(string fname)
   //{
   //   FILE *file;
   //   file = fopen(fname.c_str(), "r");

   //   if (file == NULL)
   //   {
   //      print("can not open " + fname);
   //      return;
   //   }

   //   char str[100];
   //   double tz = 0;
   //   double tvx = 0;
   //   double tvy = 0;
   //   double tvz = 0;

   //   int r = fscanf(file, "%s", str);
   //   do
   //   {
   //      r = fscanf(file, "%lf;%lf;%lf%;%lf", &tz, &tvx, &tvy, &tvz);
   //      z.push_back(tz);
   //      vx.push_back(tvx);
   //      vy.push_back(tvy);
   //      vz.push_back(tvz);
   //   } while (r != EOF);

   //   fclose(file);
   //}
   //////////////////////////////////////////////////////////////////////////
   void CalculateMeanVelocity(vector<string> &fname, string output)
   {
      print("start calculate mean velocity");
      FILE *file;
      file = fopen(fname[0].c_str(), "r");

      if (file == NULL)
      {
         print("can not open " + fname[0]);
         return;
      }

      vector<double> z;
      char str[100];
      double tz = 0;
      double tvx = 0;
      double tvy = 0;
      double tvz = 0;
      int c = 0;
      int r = fscanf(file, "%s", str);

      while (r != EOF)
      {
         r = fscanf(file, "%lf;%lf;%lf;%lf", &tz, &tvx, &tvy, &tvz);
         z.push_back(tz);
         c++;
      } 

      fclose(file);

      vector<double> avx(c);
      vector<double> avy(c);
      vector<double> avz(c);

      int ssize = fname.size();

      for (int j = 0; j < ssize; j++)
      {
         print("start read: "+fname[j]);
         FILE *file;
         file = fopen(fname[j].c_str(), "r");

         if (file == NULL)
         {
            print("can not open " + fname[j]);
            return;
         }

         vector<double> vx;
         vector<double> vy;
         vector<double> vz;

         r = fscanf(file, "%s", str);
         while (r != EOF)
         {
            r = fscanf(file, "%lf;%lf;%lf;%lf", &tz, &tvx, &tvy, &tvz);
            vx.push_back(tvx);
            vy.push_back(tvy);
            vz.push_back(tvz);
         } 

         fclose(file);

         print("end read: "+fname[j]);

         int size = (int)z.size();

#pragma omp parallel for
         for (int i = 0; i < size; i++)
         {
            avx[i] += vx[i];
            avy[i] += vy[i];
            avz[i] += vz[i];
         }
      }

#pragma omp parallel for
      for (int i = 0; i < c; i++)
      {
         avx[i] /= ssize;
         avy[i] /= ssize;
         avz[i] /= ssize;
      }

      print("start write : "+output);

      file = fopen(output.c_str(), "w");

      if (file == NULL)
      {
         print("can not open " + fname[0]);
         return;
      }

      fprintf(file, "z;<Vx>;<Vy>;<Vz>\n");

      for (int i = 0; i < c; i++)
      {
         fprintf(file, "%g;%g;%g;%g\n", z[i], avx[i], avy[i], avz[i]);
      }

      fclose(file);

      print("end write : "+output);
      print("end calculate mean velocity");
   }
   //////////////////////////////////////////////////////////////////////////
   void CalculateFluctuations(vector<string> &fname, string mv, string output)
   {
      print("start calculate fluctuations");
      FILE *file;
      file = fopen(mv.c_str(), "r");

      if (file == NULL)
      {
         print("can not open " + mv);
         return;
      }

      vector<double> z;
      vector<double> avx;
      vector<double> avy;
      vector<double> avz;

      char str[100];
      double tz = 0;
      double tvx = 0;
      double tvy = 0;
      double tvz = 0;
      int c = 0;
      int r = fscanf(file, "%s", str);

      //do
      //{
      //   r = fscanf(file, "%lf;%lf;%lf;%lf", &tz, &tvx, &tvy, &tvz);
      //   z.push_back(tz);
      //   avx.push_back(tvx);
      //   avy.push_back(tvy);
      //   avz.push_back(tvz);
      //   c++;
      //} while (r != EOF);

      while (r != EOF)
      {
         r = fscanf(file, "%lf;%lf;%lf;%lf\n", &tz, &tvx, &tvy, &tvz);
         z.push_back(tz);
         avx.push_back(tvx);
         avy.push_back(tvy);
         avz.push_back(tvz);
         c++;
      } 

      fclose(file);

      //vector<double> tfx(c,0);
      //vector<double> tfy(c,0);
      //vector<double> tfz(c,0);
      vector<double> vxx(c,0);
      vector<double> vyy(c,0);
      vector<double> vzz(c,0);
      vector<double> vxy(c,0);
      vector<double> vxz(c,0);
      vector<double> vyz(c,0);
      vector<double> vxxx(c,0);
      vector<double> vxxy(c,0);
      vector<double> vxxz(c,0);
      vector<double> vyyy(c,0);
      vector<double> vyyx(c,0);
      vector<double> vyyz(c,0);
      vector<double> vzzz(c,0);
      vector<double> vzzx(c,0);
      vector<double> vzzy(c,0);
      vector<double> vxyz(c,0);

      int ssize = fname.size();

      for (int j = 0; j < ssize; j++)
      {
         print("start read: "+fname[j]);
         FILE *file;
         file = fopen(fname[j].c_str(), "r");

         if (file == NULL)
         {
            print("can not open " + fname[j]);
            return;
         }

         vector<double> vx;
         vector<double> vy;
         vector<double> vz;

         r = fscanf(file, "%s", str);
         //do
         //{
         //   r = fscanf(file, "%lf;%lf;%lf;%lf", &tz, &tvx, &tvy, &tvz);
         //   vx.push_back(tvx);
         //   vy.push_back(tvy);
         //   vz.push_back(tvz);
         //} while (r != EOF);

         while (r != EOF)
         {
            r = fscanf(file, "%lf;%lf;%lf;%lf", &tz, &tvx, &tvy, &tvz);
            vx.push_back(tvx);
            vy.push_back(tvy);
            vz.push_back(tvz);
         }

         fclose(file);

         print("end read: "+fname[j]);

         int size = (int)z.size();

#pragma omp parallel for
         for (int i = 0; i < size; i++)
         {
            double tfx = vx[i]-avx[i];
            double tfy = vy[i]-avy[i];
            double tfz = vz[i]-avz[i];

            vxx[i] += tfx*tfx;
            vyy[i] += tfy*tfy;
            vzz[i] += tfz*tfz;
            vxy[i] += tfx*tfy;
            vxz[i] += tfx*tfz;
            vyz[i] += tfy*tfz;

            vxxx[i] += tfx*tfx*tfx;
            vxxy[i] += tfx*tfx*tfy;
            vxxz[i] += tfx*tfx*tfz;
            vyyy[i] += tfy*tfy*tfy;
            vyyx[i] += tfy*tfy*tfx;
            vyyz[i] += tfy*tfy*tfz;
            vzzz[i] += tfz*tfz*tfz;
            vzzx[i] += tfz*tfz*tfx;
            vzzy[i] += tfz*tfz*tfy;
            vxyz[i] += tfx*tfy*tfz;
         }
      }

#pragma omp parallel for
      for (int i = 0; i < c; i++)
      {
         vxx[i] /= ssize;
         vyy[i] /= ssize;
         vzz[i] /= ssize;
         vxy[i] /= ssize;
         vxz[i] /= ssize;
         vyz[i] /= ssize;

         vxxx[i] /= ssize;
         vxxy[i] /= ssize;
         vxxz[i] /= ssize;
         vyyy[i] /= ssize;
         vyyx[i] /= ssize;
         vyyz[i] /= ssize;
         vzzz[i] /= ssize;
         vzzx[i] /= ssize;
         vzzy[i] /= ssize;
         vxyz[i] /= ssize;
      }

      print("start write : "+output);

      file = fopen(output.c_str(), "w");

      if (file == NULL)
      {
         print("can not open " + fname[0]);
         return;
      }

      fprintf(file, "z;Vx;Vy;Vz;Vxx;Vyy;Vzz;Vxy;Vxz;Vyz;Vxxx;Vxxy;Vxxz;Vyyy;Vyyx;Vyyz;Vzzz;Vzzx;Vzzy;Vxyz\n");

      for (int i = 0; i < c; i++)
      {
         fprintf(file, "%g;%g;%g;%g;", z[i], avx[i], avy[i], avz[i]);
         fprintf(file, "%g;%g;%g;%g;%g;%g;", vxx[i], vyy[i], vzz[i], vxy[i], vxz[i], vyz[i]);
         fprintf(file, "%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n", vxxx[i], vxxy[i], vxxz[i], vyyy[i], vyyx[i], vyyz[i], vzzz[i], vzzx[i], vzzy[i], vxyz[i]);
      }                                 

      fclose(file);

      print("end write : "+output);
      print("end calculate fluctuations");
   }
   //////////////////////////////////////////////////////////////////////////
   void AvAllRun(char c)
   {
      using namespace std;
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();

      //string dataName = "d:/Projects/SFB880/Platte/Platte2/avs_slip84000/avs_slip84000.pvtp";
      //string output = "d:/Projects/SFB880/Platte/Platte2/avs_slip84000.csv";
      //string dataName = "d:/Projects/SFB880/Platte/PorPlatte/avs73000/avs73000.pvtp";
      //string output = "d:/Projects/SFB880/Platte/PorPlatte/avs73000/avs73000_test.csv";


      //double start = 5.2899;
      //double step = 0.078125;

      //double step = 1.0;  0.0390625;
      ////double step = 0.1;
      //double start = 10.5-step*0.5;


      vector<double> levelCoords;

      string path, path2;
      double start, step, corseDeltax;
      int minLevel, maxLevel, fstart, fstop, fstep;

      switch (c)
      {
      case 's':
         path = "/hpc3lustre/work/koskuche/SFB880/ForISM/smooth/";
         step = 0.078125;
         start = 5.25+step*0.5;
         levelCoords.push_back(55.2); //start coordinate for level 0
         levelCoords.push_back(25.2); //level 1
         levelCoords.push_back(15.2); //level 2
         levelCoords.push_back(7.75); //level 3
         levelCoords.push_back(5.25); //level 4
         minLevel = 0;
         maxLevel = 4;
         corseDeltax = 1.25;
         break;
      case 'p':
         //path = "/hpc3lustre/work/koskuche/SFB880/ForISM/porous/";
         path = "/hpc3lustre/work/koskuche/SFB880/porplate/ppv5";
         step = 0.0390625;
         //start = 10.5+step*0.5;
         start = 8.63+step*0.5;
         levelCoords.push_back(60.5); //start coordinate for level 0
         levelCoords.push_back(30.5); //level 1
         levelCoords.push_back(20.5); //level 2
         levelCoords.push_back(15.5); //level 3
         levelCoords.push_back(11.8); //level 4
         levelCoords.push_back(8.63); //level 5
         minLevel = 0;
         maxLevel = 5;
         corseDeltax = 1.25;
         fstart = 61000;
         fstop = 73000;
         //fstart = 67000;
         //fstop = 67000;
         fstep = 1000;
         break;
      }

      //for (int f = fstart; f <= fstop; f += fstep)
      //int f = 67000;
      //{
      //   path2 = path + "/pp" + toString(f);
      //   //for (int i = -10; i < 70; i += 5)
      //   for (int i = 910; i < 985; i += 5)
      //   {
      //      timer->StartTimer();
      //      print("Create averaging data for avv"+toString(i));
      //      //SurfaceAveragingAll(path+"avv"+toString(i)+"/avv"+toString(i)+".pvtu", start, step, path+"avv"+toString(i)+".csv");
      //      //SurfaceAveragingAllWithLevels(path+"avv"+toString(i)+"/avv"+toString(i)+".pvtu", minLevel, maxLevel, levelCoords, corseDeltax, path+"2/avv"+toString(i)+".csv");
      //      //SurfaceAveragingAllWithLevels(path+"/pp"+toString(i)+".pvtu", minLevel, maxLevel, levelCoords, corseDeltax, path+"/avv"+toString(i)+".csv");
      //      SurfaceAveragingVelocityWithLevels(path2+"/pp"+toString(i)+".pvtu", minLevel, maxLevel, levelCoords, corseDeltax, path2+"/avv"+toString(i)+".csv");
      //      timer->StopTimer();
      //      print("total time: " + toString(timer->GetElapsedTime()) + " s");
      //   }
      //}

      //ReadCSV("d:/Projects/SFB880/Platte/PorPlatte/TC/avv910.csv");

      double nstart = 910;
      double nstop = 985;
      double nstep = 5;

      vector<string> fstring;

      for (int i = nstart; i < nstop; i += nstep)
      {
         //for (int f = fstart; f <= fstop; f += fstep)
         //{
         //   fstring.push_back(path + "/pp" + toString(f)+"/avv"+toString(i)+".csv");
         //}
         //print("calculate mean velocity for: "+toString(i));
         //CalculateMeanVelocity(fstring, path+"/mv/mv"+toString(i)+".csv");
         //CalculateFluctuations(fstring, path+"/mv/mv"+toString(i)+".csv", path+"/cor/cor"+toString(i)+".csv");
         for (int f = fstart; f <= fstop; f += fstep)
         {
            fstring.push_back(path + "/pp" + toString(f)+"/pp"+toString(i)+".pvtu");
         }
         //print("calculate mean velocity for: "+toString(i));
         //PP::CalculateMeanVelocity(fstring, path+"/mvt/mv"+toString(i));
         //PP::CalculateFluctuations(fstring, path+"/mvt/mv"+toString(i)+".vtu", path+"/tft/tf"+toString(i));
         //SurfaceAveragingVelocityWithLevels(path+"/mvt/mv"+toString(i)+".vtu", minLevel, maxLevel, levelCoords, corseDeltax, path+"/ta/mv"+toString(i)+".csv");
         SurfaceAveragingFluctuationsWithLevels(path+"/tft/tf"+toString(i)+".vtu", minLevel, maxLevel, levelCoords, corseDeltax, path+"/ta/tf"+toString(i)+".csv");
         fstring.resize(0);
      }


      //SurfaceAveragingAll(dataName, start, step, output);
   }

}


