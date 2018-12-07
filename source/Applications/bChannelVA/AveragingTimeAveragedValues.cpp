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

namespace ATAV
{

   //////////////////////////////////////////////////////////////////////////
   void SurfaceAveragingAllWithLevels(std::string dataName, vector<int>& levels, vector<double>& levelCoords, double corseDeltax, std::string output)
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

      fprintf(file, "z;<Vx>;<Vy>;<Vz>;<Vxx>;<Vyy>;<Vzz>;<Vxy>;<Vxz>;<Vyz>;<Vxxx>;<Vxxy>;<Vxxz>;<Vyyy>;<Vyyx>;<Vyyz>;<Vzzz>;<Vzzx>;<Vzzy>;<Vxyz>\n");

      int size = (int)levels.size();

      for (int i = 0; i < size; i++)
      {
         int level = levels[i];
         double step = corseDeltax/(double)(1<<level);
         double start = levelCoords[i]+step*0.5;
         double stop = levelCoords[i + 1] - step*0.5;

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
            vtkSmartPointer<vtkDataArray> vxArray = planeCut->GetOutput()->GetPointData()->GetArray("taVx");
            vtkSmartPointer<vtkDataArray> vyArray = planeCut->GetOutput()->GetPointData()->GetArray("taVy");
            vtkSmartPointer<vtkDataArray> vzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVz");

            vtkSmartPointer<vtkDataArray> vxxArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxx");
            vtkSmartPointer<vtkDataArray> vyyArray = planeCut->GetOutput()->GetPointData()->GetArray("taVyy");
            vtkSmartPointer<vtkDataArray> vzzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVzz");
            vtkSmartPointer<vtkDataArray> vxyArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxy");
            vtkSmartPointer<vtkDataArray> vxzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxz");
            vtkSmartPointer<vtkDataArray> vyzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVyz");

            vtkSmartPointer<vtkDataArray> vxxxArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxxx");
            vtkSmartPointer<vtkDataArray> vxxyArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxxy");
            vtkSmartPointer<vtkDataArray> vxxzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxxz");
            vtkSmartPointer<vtkDataArray> vyyyArray = planeCut->GetOutput()->GetPointData()->GetArray("taVyyy");
            vtkSmartPointer<vtkDataArray> vyyxArray = planeCut->GetOutput()->GetPointData()->GetArray("taVyyx");
            vtkSmartPointer<vtkDataArray> vyyzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVyyz");
            vtkSmartPointer<vtkDataArray> vzzzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVzzz");
            vtkSmartPointer<vtkDataArray> vzzxArray = planeCut->GetOutput()->GetPointData()->GetArray("taVzzx");
            vtkSmartPointer<vtkDataArray> vzzyArray = planeCut->GetOutput()->GetPointData()->GetArray("taVzzy");
            vtkSmartPointer<vtkDataArray> vxyzArray = planeCut->GetOutput()->GetPointData()->GetArray("taVxyz");

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
               vxAvSumm += vxArray->GetTuple1(i);
               vyAvSumm += vyArray->GetTuple1(i);
               vzAvSumm += vzArray->GetTuple1(i);

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
               vxAv = vxAvSumm/numPoints;
               vyAv = vyAvSumm/numPoints;
               vzAv = vzAvSumm/numPoints;

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

            fprintf(file, "%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n", j, vxAv, vyAv, vzAv, vxxAv, vyyAv, vzzAv, vxyAv, vxzAv, vyzAv,
               vxxxAv, vxxyAv, vxxzAv, vyyyAv, vyyxAv, vyyzAv, vzzzAv, vzzxAv, vzzyAv, vxyzAv);

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
      vector<int> levels;

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
         path = "/hpc3lustre/work/koskuche/SFB880/porplate/avv/96000";
         
         levelCoords.push_back(8.63); //level 5
         levelCoords.push_back(11.8); //level 4
         levelCoords.push_back(15.5); //level 3
         levelCoords.push_back(20.5); //level 2
         levelCoords.push_back(30.5); //level 1
         levelCoords.push_back(60.5); //start coordinate for level 0
         levelCoords.push_back(300.0);// z max

         levels.push_back(5);
         levels.push_back(4);
         levels.push_back(3);
         levels.push_back(2);
         levels.push_back(1);
         levels.push_back(0);

         corseDeltax = 1.25;
         break;
      case 'c':
         path = "/hpc3lustre/work/koskuche/SFB880/pChannel/pChannelHLRN/va/120000";
         levelCoords.push_back(0); //start coordinate for level 2
         levelCoords.push_back(0.004); //level 1
         levelCoords.push_back(0.006); //level 0
         levelCoords.push_back(0.016); //level 1
         levelCoords.push_back(0.019); //level 2
         levelCoords.push_back(0.02); //z max
         
         levels.push_back(2);
         levels.push_back(1);
         levels.push_back(0);
         levels.push_back(1);
         levels.push_back(2);

         corseDeltax = 8e-05;
         break;
      }

      int i=35;
      //SurfaceAveragingAllWithLevels(path+"avv"+toString(i)+"/avv"+toString(i)+".pvtu", levels, levelCoords, corseDeltax, path+"avv"+toString(i)+".csv");
      //SurfaceAveragingAllWithLevels(path + "/va.vtu", levels, levelCoords, corseDeltax, path + "/avv1.csv");
      SurfaceAveragingAllWithLevels(path + "/avv35.vtu", levels, levelCoords, corseDeltax, path + "/avv96.csv");
   }

}


