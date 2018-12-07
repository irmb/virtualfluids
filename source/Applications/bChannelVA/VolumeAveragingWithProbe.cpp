#include "Postprocessing.h"
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkTimerLog.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkThreshold.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkFileOutputWindow.h>
#include "MemoryUtil.h"
#include <algorithm>
#include <numeric>
using namespace std;

namespace VAP
{
   double deltax;
   double dx3;
   double l;
   double lQuadrat;
   double lNorm;

   int geo_nx1, geo_nx2, geo_nx3;
   int geo_extent[6];
   double geo_origin[3];
   double geo_spacing[3];

   double bbox[6];

   vector<int> geoMatrix;
   vector<int> lvlMatrix;
   vector<double> vxMatrix;
   vector<double> vyMatrix;
   vector<double> vzMatrix;
   vector<double> prMatrix;
   vector<double> vaVxMatrix;
   vector<double> vaVyMatrix;
   vector<double> vaVzMatrix;
   vector<double> vaPrMatrix;

   //////////////////////////////////////////////////////////////////////////
   inline int index(int x1, int x2, int x3)
   {
      return geo_nx1 * (geo_nx2 * x3 + x2) + x1;
   }
   //////////////////////////////////////////////////////////////////////////
   inline int index(int x1, int x2, int x3, int nx1, int nx2)
   {
      return nx1 * (nx2 * x3 + x2) + x1;
   }

   //////////////////////////////////////////////////////////////////////////
   void getNodeIndexes(double x[3], int ix[3], double origin[3], double deltax)
   {
      ix[0] = cint((x[0]-origin[0])/deltax);
      ix[1] = cint((x[1]-origin[1])/deltax);
      ix[2] = cint((x[2]-origin[2])/deltax);
   }
   //////////////////////////////////////////////////////////////////////////
   void getNodeCoordinates(double x[3], int ix[3], double origin[3], double deltax)
   {
      x[0] = origin[0] + ix[0]*deltax;
      x[1] = origin[1] + ix[1]*deltax;
      x[2] = origin[2] + ix[2]*deltax;
   }
   //////////////////////////////////////////////////////////////////////////
   double trilinearInterpolation(double x[3], double y[3], double z[3], double val[8])
   {
      double c = (x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1]);

      return val[0]/c*(x[2]-x[0])*(y[2]-y[0])*(z[2]-z[0])+
         val[1]/c*(x[2]-x[0])*(y[2]-y[0])*(z[0]-z[1])+
         val[2]/c*(x[2]-x[0])*(y[0]-y[1])*(z[2]-z[0])+
         val[3]/c*(x[2]-x[0])*(y[0]-y[1])*(z[0]-z[1])+
         val[4]/c*(x[0]-x[1])*(y[2]-y[0])*(z[2]-z[0])+
         val[5]/c*(x[0]-x[1])*(y[2]-y[0])*(z[0]-z[1])+
         val[6]/c*(x[0]-x[1])*(y[0]-y[1])*(z[2]-z[0])+
         val[7]/c*(x[0]-x[1])*(y[0]-y[1])*(z[0]-z[1]);
   }
   //////////////////////////////////////////////////////////////////////////
   void initBoundingBox(string blocksName)
   {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + blocksName + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetBlocks(ReadDataSet(blocksName.c_str()));
      print("read data set from " + blocksName + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      dataSetBlocks->GetBounds(bbox);
      print("Bounding Box: NX1 x NX2 x NX3 = " + toString(bbox[1]) + " x " + toString(bbox[3]) + " x " + toString(bbox[5]));

      //geo_origin[0] = bbox[0]-0.5*deltax;
      //geo_origin[1] = bbox[2]-0.5*deltax;
      //geo_origin[2] = bbox[4]-0.5*deltax;

      //geo_nx1 = cint((bbox[1] - bbox[0]) / deltax+2);
      //geo_nx2 = cint((bbox[3] - bbox[2]) / deltax+2);
      //geo_nx3 = cint((bbox[5] - bbox[4]) / deltax+2);

      geo_origin[0] = bbox[0]+0.5*deltax;
      geo_origin[1] = bbox[2]+0.5*deltax;
      geo_origin[2] = bbox[4]+0.5*deltax;

      geo_nx1 = cint((bbox[1] - bbox[0]) / deltax-2);
      geo_nx2 = cint((bbox[3] - bbox[2]) / deltax-2);
      geo_nx3 = cint((bbox[5] - bbox[4]) / deltax-2);

      print("geo NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      geo_extent[0] = 0;
      geo_extent[2] = 0;
      geo_extent[4] = 0;
      geo_extent[1] = geo_nx1-1;
      geo_extent[3] = geo_nx2-1;
      geo_extent[5] = geo_nx3-1;

      geo_spacing[0] = deltax;
      geo_spacing[1] = deltax;
      geo_spacing[2] = deltax;

      int size = geo_nx1*geo_nx2*geo_nx3;
      lvlMatrix.resize(size, 0);
      //geoMatrix.resize(size, 0);

      vtkSmartPointer<vtkDataArray> lArray = dataSetBlocks->GetCellData()->GetArray("level");

      double range[2];
      lArray->GetRange(range);

      for (int level = (int)range[0]; level <=(int)range[1]; level++)
      {
         print("Perform the blocks level " + toString(level) + " : start");
         level_interp_timer->StartTimer();

         vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
         thrFilter->SetInputData(dataSetBlocks);
         thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "level");
         thrFilter->ThresholdBetween(level, level);
         thrFilter->Update();
         vtkSmartPointer<vtkUnstructuredGrid> lugrid = thrFilter->GetOutput();

         vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
         vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

         int numberOfCells = lugrid->GetNumberOfCells();

         double x[3];
         double xMin[3];
         double xMax[3];
         int ixMin[3];
         int ixMax[3];
         vtkIdType idc = 0;

         for (int i = 0; i < numberOfCells; i++)
         {
            vtkIdType npts = 8;
            vtkIdType* pts;
            lugrid->GetCellPoints(i, npts, pts);
            vector <double> x1;
            vector <double> x2;
            vector <double> x3;
            for (int p = 0; p < 8; p++)
            {
               lugrid->GetPoint(pts[p], x);
               x1.push_back(x[0]);
               x2.push_back(x[1]);
               x3.push_back(x[2]);
            }
            xMin[0] = *min_element(x1.begin(), x1.end());
            xMin[1] = *min_element(x2.begin(), x2.end());
            xMin[2] = *min_element(x3.begin(), x3.end());

            xMax[0] = *max_element(x1.begin(), x1.end());
            xMax[1] = *max_element(x2.begin(), x2.end());
            xMax[2] = *max_element(x3.begin(), x3.end());

            getNodeIndexes(xMin, ixMin, geo_origin, deltax);
            getNodeIndexes(xMax, ixMax, geo_origin, deltax);

            for (int k = ixMin[2]; k <= ixMax[2]; k++)
            {
               for (int j = ixMin[1]; j <= ixMax[1]; j++)
               {
                  for (int i = ixMin[0]; i <= ixMax[0]; i++)
                  {
                     if (i >= 0 && j >= 0 && k >= 0 && i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                     {
                        lvlMatrix[index(i, j, k)] = level;
                     }
                  }
               }
            }
         }
         print("Perform the blocks level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }

   }
   //////////////////////////////////////////////////////////////////////////
   void createGeoMatrix(string dataNameG)
   {
      print("createGeoMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      //print("read data set from " + dataName + ": start");
      //timer->StartTimer();
      //vtkSmartPointer<vtkDataSet> dataSet(ReadDataSet(dataName.c_str()));
      //print("read data set from " + dataName + ": end");
      //timer->StopTimer();
      //print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      double bounds[6];
      dataSetGeo->GetBounds(bounds);

      double origin[3];
      origin[0] = bounds[0] + deltax / 2.0;
      origin[1] = bounds[2] + deltax / 2.0;
      origin[2] = bounds[4] + deltax / 2.0;

      geo_nx1 = int((bounds[1] - bounds[0]) / deltax);
      geo_nx2 = int((bounds[3] - bounds[2]) / deltax);
      geo_nx3 = int((bounds[5] - bounds[4]) / deltax);

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));


      int size = geo_nx1*geo_nx2*geo_nx3;
      geoMatrix.resize(size, 0);

      geo_extent[0] = 0;
      geo_extent[2] = 0;
      geo_extent[4] = 0;
      geo_extent[1] = geo_nx1 - 1;
      geo_extent[3] = geo_nx2 - 1;
      geo_extent[5] = geo_nx3 - 1;

      geo_origin[0] = origin[0];
      geo_origin[1] = origin[1];
      geo_origin[2] = origin[2];

      geo_spacing[0] = deltax;
      geo_spacing[1] = deltax;
      geo_spacing[2] = deltax;

      int blocknx = 100;
      double blockln = (double)(blocknx)*deltax;
      int i = 0;

      int partx1 = (int)ceil((double)geo_nx1 / (double)(blocknx));
      int partx2 = (int)ceil((double)geo_nx2 / (double)(blocknx));
      int partx3 = (int)ceil((double)geo_nx3 / (double)(blocknx));

      int maxX1 = blocknx*(int)partx1;
      int maxX2 = blocknx*(int)partx2;
      int maxX3 = blocknx*(int)partx3;

      vtkSmartPointer<vtkImageData> vmatrix = vtkImageData::New();
      int extent[6];
      extent[0] = 0;
      extent[2] = 0;
      extent[4] = 0;
      extent[1] = blocknx - 1;
      extent[3] = blocknx - 1;
      extent[5] = blocknx - 1;
      vmatrix->SetExtent(extent);
      vmatrix->SetSpacing(deltax, deltax, deltax);

      // Perform the interpolation
      vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
      probeFilter->SetSourceData(dataSetGeo);
      // Interpolate 'Source' at these points
#if VTK_MAJOR_VERSION <= 5
      probeFilter->SetInput(image);
#else
      probeFilter->SetInputData(vmatrix);
#endif


      for (int x3 = 0; x3 < partx3; x3++)
      {
         for (int x2 = 0; x2 < partx2; x2++)
         {
            for (int x1 = 0; x1 < partx1; x1++)
            {
               vmatrix->SetOrigin(origin[0] + x1*blockln, origin[1] + x2*blockln, origin[2] + x3*blockln);
               print("start interpolation");
               timer->StartTimer();

               probeFilter->Update();
               timer->StopTimer();
               print("end interpolation");
               print("interpolation time: " + toString(timer->GetElapsedTime()) + " s");

               int id = 0;

               for (int k = x3*blocknx; k < x3*blocknx + blocknx; k++)
               {
                  for (int j = x2*blocknx; j < x2*blocknx + blocknx; j++)
                  {
                     for (int i = x1*blocknx; i < x1*blocknx + blocknx; i++)
                     {
                        int g = (int)probeFilter->GetOutput()->GetPointData()->GetArray("Geometry")->GetTuple1(id);
                        if (i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                        {
                           if (g == 1)
                              geoMatrix[index(i, j, k)] = 0;
                           else
                              geoMatrix[index(i, j, k)] = 1;
                        }
                        id++;
                     }
                  }
               }

               print("i=" + toString(i));
               i++;
            }
         }
      }

      print("createGeoMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void readGeoMatrix(string dataNameG)
   {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();

      print("readGeoMatrix:start");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      vtkSmartPointer<vtkImageData> image = vtkImageData::SafeDownCast(dataSetGeo);

      image->GetExtent(geo_extent);
      image->GetOrigin(geo_origin);
      image->GetSpacing(geo_spacing);

      geo_nx1 = geo_extent[1] + 1;
      geo_nx2 = geo_extent[3] + 1;
      geo_nx3 = geo_extent[5] + 1;

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      int size = geo_nx1*geo_nx2*geo_nx3;
      geoMatrix.resize(size, 0);
      lvlMatrix.resize(size, 0);

      vtkSmartPointer<vtkDataArray> geoArray = dataSetGeo->GetPointData()->GetArray("geo");
      vtkSmartPointer<vtkDataArray> lvlArray = dataSetGeo->GetPointData()->GetArray("level");

      int numberOfPoints = dataSetGeo->GetNumberOfPoints();

      for (int i = 0; i < numberOfPoints; i++)
      {
         geoMatrix[i] = (int)geoArray->GetTuple1(i);
         lvlMatrix[i] = (int)lvlArray->GetTuple1(i);
      }

      print("readGeoMatrix:end");
   }
   //////////////////////////////////////////////////////////////////////////
   void createMQMatrix(string dataNameMQ, std::string dataNameG)
   {
      print("createMQMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetMQ(ReadDataSet(dataNameMQ.c_str()));
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");
      vtkSmartPointer<vtkDataArray> vxArray = dataSetMQ->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSetMQ->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSetMQ->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Press");

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      int blocknx = 100;
      double blockln = (double)(blocknx)*deltax;
      int count = 0;

      int partx1 = (int)ceil((double)geo_nx1 / (double)(blocknx));
      int partx2 = (int)ceil((double)geo_nx2 / (double)(blocknx));
      int partx3 = (int)ceil((double)geo_nx3 / (double)(blocknx));

      int maxX1 = blocknx*(int)partx1;
      int maxX2 = blocknx*(int)partx2;
      int maxX3 = blocknx*(int)partx3;

      vtkSmartPointer<vtkImageData> vmatrix = vtkImageData::New();
      int extent[6];
      extent[0] = 0;
      extent[2] = 0;
      extent[4] = 0;
      extent[1] = blocknx - 1;
      extent[3] = blocknx - 1;
      extent[5] = blocknx - 1;
      vmatrix->SetExtent(extent);
      vmatrix->SetSpacing(deltax, deltax, deltax);


      // Perform the interpolation
      vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();

      double range[2];
      lArray->GetRange(range);

      for (int level = (int)range[0]; level <= (int)range[1]; level++)
      {
         print("Perform the interpolation level " + toString(level) + " : start");
         level_interp_timer->StartTimer();

         print("actual memory usage befor interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

         print("Create the level grid " + toString(level) + " : start");
         level_grid_timer->StartTimer();

         vtkSmartPointer<vtkDoubleArray> vxArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vxArrayNew->SetNumberOfComponents(1);
         vxArrayNew->SetName("Vx");

         vtkSmartPointer<vtkDoubleArray> vyArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vyArrayNew->SetNumberOfComponents(1);
         vyArrayNew->SetName("Vy");

         vtkSmartPointer<vtkDoubleArray> vzArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vzArrayNew->SetNumberOfComponents(1);
         vzArrayNew->SetName("Vz");

         vtkSmartPointer<vtkDoubleArray> prArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         prArrayNew->SetNumberOfComponents(1);
         prArrayNew->SetName("Press");

         vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

         vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
         vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

         int numberOfCells = dataSetMQ->GetNumberOfCells();

         double x[3];
         vtkIdType idc = 0;

         for (int i = 0; i < numberOfCells; i++)
         {
            vtkSmartPointer<vtkIdList> plist = vtkSmartPointer<vtkIdList>::New();
            dataSetMQ->GetCellPoints(i, plist);
            dataSetMQ->GetPoint(plist->GetId(0), x);
            vtkIdType id = dataSetGeo->FindPoint(x);
            if (!((int)id < 0))
            {
               int l = (int)lArray->GetTuple1(id);
               if (l == level)
               {
                  vtkSmartPointer<vtkIdList> plist2 = vtkSmartPointer<vtkIdList>::New();
                  for (int p = 0; p < plist->GetNumberOfIds(); p++)
                  {
                     dataSetMQ->GetPoint(plist->GetId(p), x);
                     points->InsertNextPoint(x);
                     vxArrayNew->InsertNextTuple1(vxArray->GetTuple1(plist->GetId(p)));
                     vyArrayNew->InsertNextTuple1(vyArray->GetTuple1(plist->GetId(p)));
                     vzArrayNew->InsertNextTuple1(vzArray->GetTuple1(plist->GetId(p)));
                     prArrayNew->InsertNextTuple1(prArray->GetTuple1(plist->GetId(p)));
                     plist2->InsertNextId(idc);
                     idc++;
                  }
                  cells->InsertNextCell(plist2);
               }
            }
         }

         ugrid->SetPoints(points);
         ugrid->SetCells(VTK_VOXEL, cells);
         ugrid->GetPointData()->AddArray(vxArrayNew);
         ugrid->GetPointData()->AddArray(vyArrayNew);
         ugrid->GetPointData()->AddArray(vzArrayNew);
         ugrid->GetPointData()->AddArray(prArrayNew);

         print("Create the level grid " + toString(level) + " : end");
         level_grid_timer->StopTimer();
         print("creating time: " + toString(level_grid_timer->GetElapsedTime()) + " s");

         probeFilter->SetSourceData(ugrid);
         //probeFilter->SetSourceData(dataSetMQ);
         // Interpolate 'Source' at these points
#if VTK_MAJOR_VERSION <= 5
         probeFilter->SetInput(vmatrix);
#else
         probeFilter->SetInputData(vmatrix);
#endif

         int count = 0;

         for (int x3 = 0; x3 < partx3; x3++)
         {
            for (int x2 = 0; x2 < partx2; x2++)
            {
               for (int x1 = 0; x1 < partx1; x1++)
               {
                  vmatrix->SetOrigin(geo_origin[0] + x1*blockln, geo_origin[1] + x2*blockln, geo_origin[2] + x3*blockln);
                  print("start interpolation");
                  timer->StartTimer();

                  probeFilter->Update();
                  timer->StopTimer();
                  print("end interpolation");
                  print("interpolation time: " + toString(timer->GetElapsedTime()) + " s");

                  int id = 0;

                  vtkSmartPointer<vtkDataArray> vxArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vx");
                  vtkSmartPointer<vtkDataArray> vyArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vy");
                  vtkSmartPointer<vtkDataArray> vzArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vz");
                  vtkSmartPointer<vtkDataArray> prArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Press");

                  for (int k = x3*blocknx; k < x3*blocknx + blocknx; k++)
                  {
                     for (int j = x2*blocknx; j < x2*blocknx + blocknx; j++)
                     {
                        for (int i = x1*blocknx; i < x1*blocknx + blocknx; i++)
                        {
                           if (i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                           {
                              double vx = vxArrayP->GetTuple1(id);
                              if (vx != 0)
                              {
                                 vxMatrix[index(i, j, k)] = vx;
                                 vyMatrix[index(i, j, k)] = vyArrayP->GetTuple1(id);
                                 vzMatrix[index(i, j, k)] = vzArrayP->GetTuple1(id);
                                 prMatrix[index(i, j, k)] = prArrayP->GetTuple1(id);
                              }
                           }
                           id++;
                        }
                     }
                  }

                  print("count=" + toString(count));
                  count++;
               }
            }
         }
         print("Perform the interpolation level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }

      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   void createMQMatrix2(string dataNameMQ, std::string dataNameG, double deltax_corse)
   {
      print("createMQMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetMQ(ReadDataSet(dataNameMQ.c_str()));
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");
      vtkSmartPointer<vtkDataArray> vxArray = dataSetMQ->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSetMQ->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSetMQ->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Press");

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      int blocknx = 100;
      double blockln = (double)(blocknx)*deltax;
      int count = 0;

      int partx1 = (int)ceil((double)geo_nx1 / (double)(blocknx));
      int partx2 = (int)ceil((double)geo_nx2 / (double)(blocknx));
      int partx3 = (int)ceil((double)geo_nx3 / (double)(blocknx));

      int maxX1 = blocknx*(int)partx1;
      int maxX2 = blocknx*(int)partx2;
      int maxX3 = blocknx*(int)partx3;

      vtkSmartPointer<vtkImageData> vmatrix = vtkImageData::New();
      int extent[6];
      extent[0] = 0;
      extent[2] = 0;
      extent[4] = 0;
      extent[1] = blocknx - 1;
      extent[3] = blocknx - 1;
      extent[5] = blocknx - 1;
      vmatrix->SetExtent(extent);
      vmatrix->SetSpacing(deltax, deltax, deltax);


      double range[2];
      lArray->GetRange(range);

      for (int level = (int)range[0]; level <=(int)range[1]; level++)
      {
         print("Perform the interpolation level " + toString(level) + " : start");
         level_interp_timer->StartTimer();

         print("actual memory usage befor interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

         print("Create the level grid " + toString(level) + " : start");
         level_grid_timer->StartTimer();

         vtkSmartPointer<vtkDoubleArray> vxArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vxArrayNew->SetNumberOfComponents(1);
         vxArrayNew->SetName("Vx");

         vtkSmartPointer<vtkDoubleArray> vyArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vyArrayNew->SetNumberOfComponents(1);
         vyArrayNew->SetName("Vy");

         vtkSmartPointer<vtkDoubleArray> vzArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vzArrayNew->SetNumberOfComponents(1);
         vzArrayNew->SetName("Vz");

         vtkSmartPointer<vtkDoubleArray> prArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         prArrayNew->SetNumberOfComponents(1);
         prArrayNew->SetName("Press");

         double bounds[6];
         vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
         thrFilter->SetInputData(dataSetGeo);
         thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Level");
         thrFilter->ThresholdBetween(level, level);
         thrFilter->Update();
         vtkSmartPointer<vtkUnstructuredGrid> lugrid = thrFilter->GetOutput();
         lugrid->GetBounds(bounds);

         vtkSmartPointer<vtkImageData> limage = vtkImageData::New();
         double delta = deltax_corse/(double)(1<<level);

         double liorigin[3];
         liorigin[0] = bounds[0];
         liorigin[1] = bounds[2];
         liorigin[2] = bounds[4];

         int nx1 = cint((bounds[1] - bounds[0]) / delta)+1;
         int nx2 = cint((bounds[3] - bounds[2]) / delta)+1;
         int nx3 = cint((bounds[5] - bounds[4]) / delta)+1;

         print("level matrix NX1 x NX2 x NX3 = " + toString(nx1) + " x " + toString(nx2) + " x " + toString(nx3));

         int liextent[6];
         liextent[0] = 0;
         liextent[2] = 0;
         liextent[4] = 0;
         liextent[1] = nx1-1;
         liextent[3] = nx2-1;
         liextent[5] = nx3-1;

         double lispacing[3];
         lispacing[0] = delta;
         lispacing[1] = delta;
         lispacing[2] = delta;

         limage->SetExtent(liextent);
         limage->SetOrigin(liorigin[0], liorigin[1], liorigin[2]);
         limage->SetSpacing(delta, delta, delta);

         int numberOfPoints = lugrid->GetNumberOfPoints();

         int msize = limage->GetNumberOfPoints();
         vector<double> lvxMatrix(msize, 0);
         vector<double> lvyMatrix(msize, 0);
         vector<double> lvzMatrix(msize, 0);
         vector<double> lprMatrix(msize, 0);

         double x[3];
         int   ix[3];

         for (int i = 0; i < numberOfPoints; i++)
         {
            lugrid->GetPoint(i, x);
            vtkIdType id = dataSetMQ->FindPoint(x);
            if (!((int)id < 0))
            {
               getNodeIndexes(x, ix, liorigin, delta);
               lvxMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = vxArray->GetTuple1(id);
               lvyMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = vyArray->GetTuple1(id);
               lvzMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = vzArray->GetTuple1(id);
               lprMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = prArray->GetTuple1(id);
            }
         }

         for (int x3 = 0; x3 < nx3; x3++)
            for (int x2 = 0; x2 < nx2; x2++)
               for (int x1 = 0; x1 < nx1; x1++)
               {
                  vxArrayNew->InsertNextValue(lvxMatrix[index(x1, x2, x3, nx1, nx2)]);
                  vyArrayNew->InsertNextValue(lvyMatrix[index(x1, x2, x3, nx1, nx2)]);
                  vzArrayNew->InsertNextValue(lvzMatrix[index(x1, x2, x3, nx1, nx2)]);
                  prArrayNew->InsertNextValue(lprMatrix[index(x1, x2, x3, nx1, nx2)]);
               }


         limage->GetPointData()->AddArray(vxArrayNew);
         limage->GetPointData()->AddArray(vyArrayNew);
         limage->GetPointData()->AddArray(vzArrayNew);
         limage->GetPointData()->AddArray(prArrayNew);

         vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
         writer->SetInputData(limage);
         writer->SetFileName("d:/temp/sfb1/limage.vti");
         writer->SetDataModeToAscii();
         writer->Update();

         print("Create the level grid " + toString(level) + " : end");
         level_grid_timer->StopTimer();
         print("creating time: " + toString(level_grid_timer->GetElapsedTime()) + " s");

         // Perform the interpolation
         vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
         probeFilter->SetSourceData(limage);

         // Interpolate 'Source' at these points
#if VTK_MAJOR_VERSION <= 5
         probeFilter->SetInput(vmatrix);
#else
         probeFilter->SetInputData(vmatrix);
#endif

         int count = 0;

         for (int x3 = 0; x3 < partx3; x3++)
         {
            for (int x2 = 0; x2 < partx2; x2++)
            {
               for (int x1 = 0; x1 < partx1; x1++)
               {
                  vmatrix->SetOrigin(geo_origin[0] + x1*blockln, geo_origin[1] + x2*blockln, geo_origin[2] + x3*blockln);
                  print("x1="+toString(x1)+" x2="+toString(x2)+" x3="+toString(x3));
                  print("origin: " + toString(geo_origin[0] + x1*blockln)+toString(geo_origin[1] + x2*blockln)+toString(geo_origin[2] + x3*blockln));
                  print("start interpolation");
                  timer->StartTimer();

                  probeFilter->Update();
                  timer->StopTimer();
                  print("end interpolation");
                  print("interpolation time: " + toString(timer->GetElapsedTime()) + " s");

                  int id = 0;

                  vtkSmartPointer<vtkDataArray> vxArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vx");
                  vtkSmartPointer<vtkDataArray> vyArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vy");
                  vtkSmartPointer<vtkDataArray> vzArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vz");
                  vtkSmartPointer<vtkDataArray> prArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Press");

                  for (int k = x3*blocknx; k < x3*blocknx + blocknx; k++)
                  {
                     for (int j = x2*blocknx; j < x2*blocknx + blocknx; j++)
                     {
                        for (int i = x1*blocknx; i < x1*blocknx + blocknx; i++)
                        {
                           if (i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                           {
                              double vx = vxArrayP->GetTuple1(id);
                              if (vx != 0)
                              {
                                 vxMatrix[index(i, j, k)] = vx;
                                 vyMatrix[index(i, j, k)] = vyArrayP->GetTuple1(id);
                                 vzMatrix[index(i, j, k)] = vzArrayP->GetTuple1(id);
                                 prMatrix[index(i, j, k)] = prArrayP->GetTuple1(id);
                              }
                           }
                           id++;
                        }
                     }
                  }

                  print("count=" + toString(count));
                  count++;
               }
            }
         }
         print("Perform the interpolation level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }

      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void createGeoMatrix3(std::string dataNameG, double deltax_corse)
   {
      print("createGeoMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");
      vtkSmartPointer<vtkDataArray> gArray = dataSetGeo->GetPointData()->GetArray("Geometry");

      //geo_origin[0] = bbox[0]+0.5*deltax;
      //geo_origin[1] = bbox[2]+0.5*deltax;
      //geo_origin[2] = bbox[4]+0.5*deltax;

      //geo_nx1 = cint((bbox[1] - bbox[0]) / deltax);
      //geo_nx2 = cint((bbox[3] - bbox[2]) / deltax);
      //geo_nx3 = cint((bbox[5] - bbox[4]) / deltax);

      //print("level matrix NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      //geo_extent[0] = 0;
      //geo_extent[2] = 0;
      //geo_extent[4] = 0;
      //geo_extent[1] = geo_nx1-1;
      //geo_extent[3] = geo_nx2-1;
      //geo_extent[5] = geo_nx3-1;

      //geo_spacing[0] = deltax;
      //geo_spacing[1] = deltax;
      //geo_spacing[2] = deltax;

      int size = geo_nx1*geo_nx2*geo_nx3;
      geoMatrix.resize(size, 1);

      vtkSmartPointer<vtkImageData> gimage = vtkSmartPointer<vtkImageData>::New();
      gimage->SetOrigin(geo_origin);
      gimage->SetExtent(geo_extent);
      gimage->SetSpacing(geo_spacing);


      double range[2];
      lArray->GetRange(range);

      for (int level = (int)range[0]; level <=(int)range[1]; level++)
      {
         print("Perform the interpolation level " + toString(level) + " : start");
         level_interp_timer->StartTimer();

         print("actual memory usage befor interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

         print("Create the level grid " + toString(level) + " : start");
         level_grid_timer->StartTimer();

         vtkSmartPointer<vtkDoubleArray> gArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         gArrayNew->SetNumberOfComponents(1);
         gArrayNew->SetName("Geometry");

         double bounds[6];
         vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
         thrFilter->SetInputData(dataSetGeo);
         thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Level");
         thrFilter->ThresholdBetween(level, level);
         thrFilter->Update();
         vtkSmartPointer<vtkUnstructuredGrid> lugrid = thrFilter->GetOutput();
         lugrid->GetBounds(bounds);

         vtkSmartPointer<vtkImageData> limage = vtkSmartPointer<vtkImageData>::New();
         double delta = deltax_corse/(double)(1<<level);

         double liorigin[3];
         liorigin[0] = bbox[0]-0.5*delta;
         liorigin[1] = bbox[2]-0.5*delta;
         liorigin[2] = bbox[4]-0.5*delta;

         int nx1 = cint((bbox[1] - bbox[0]) / deltax+2);
         int nx2 = cint((bbox[3] - bbox[2]) / deltax+2);
         int nx3 = cint((bbox[5] - bbox[4]) / deltax+2);

         print("level matrix NX1 x NX2 x NX3 = " + toString(nx1) + " x " + toString(nx2) + " x " + toString(nx3));

         int liextent[6];
         liextent[0] = 0;
         liextent[2] = 0;
         liextent[4] = 0;
         liextent[1] = nx1-1;
         liextent[3] = nx2-1;
         liextent[5] = nx3-1;

         double lispacing[3];
         lispacing[0] = delta;
         lispacing[1] = delta;
         lispacing[2] = delta;

         limage->SetExtent(liextent);
         limage->SetOrigin(liorigin[0], liorigin[1], liorigin[2]);
         limage->SetSpacing(delta, delta, delta);

         int numberOfPoints = lugrid->GetNumberOfPoints();

         int msize = limage->GetNumberOfPoints();
         vector<double> lgeoMatrix(msize, 0);

         double x[3];
         int   ix[3];

         for (int i = 0; i < numberOfPoints; i++)
         {
            lugrid->GetPoint(i, x);
            vtkIdType id = dataSetGeo->FindPoint(x);
            if (!((int)id < 0))
            {
               getNodeIndexes(x, ix, liorigin, delta);
               //print("ix0="+toString(ix[0])+" ix1="+ toString(ix[1])+" ix2="+toString(ix[2]));
               if (ix[0] >= 0 && ix[1] >= 0 && ix[2] >= 0 && ix[0] < nx1 && ix[1] < nx2 && ix[2] < nx3)
               {
                  lgeoMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = gArray->GetTuple1(id);

                  //lvlMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = level;
               }
            }
         }


         for (int x3 = 0; x3 < nx3; x3++)
            for (int x2 = 0; x2 < nx2; x2++)
               for (int x1 = 0; x1 < nx1; x1++)
               {
                  gArrayNew->InsertNextValue(lgeoMatrix[index(x1, x2, x3, nx1, nx2)]);
               }


         limage->GetPointData()->AddArray(gArrayNew);

         //vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
         //writer->SetInputData(limage);
         //writer->SetFileName("d:/temp/sfb1/limage.vti");
         //writer->SetDataModeToAscii();
         //writer->Update();

         print("Create the level grid " + toString(level) + " : end");
         level_grid_timer->StopTimer();
         print("creating time: " + toString(level_grid_timer->GetElapsedTime()) + " s");

         // Perform the interpolation
         vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
         probeFilter->SetSourceData(limage);

         // Interpolate 'Source' at these points
#if VTK_MAJOR_VERSION <= 5
         probeFilter->SetInput(gimage);
#else
         probeFilter->SetInputData(gimage);
#endif

         print("start interpolation");
         timer->StartTimer();

         probeFilter->Update();
         timer->StopTimer();
         print("end interpolation");
         print("interpolation time: " + toString(timer->GetElapsedTime()) + " s");

         int id = 0;

         vtkSmartPointer<vtkDataArray> geoArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Geometry");

         for (int k = 0; k < geo_nx3; k++)
         {
            for (int j = 0; j < geo_nx2; j++)
            {
               for (int i = 0; i < geo_nx1; i++)
               {
                  if (i < geo_nx1 && j < geo_nx2 && k < geo_nx3 && geoMatrix[index(i, j, k)] == 1)
                  {
                     int g = (int)geoArrayP->GetTuple1(id);
                     if (g == 1)
                        geoMatrix[index(i, j, k)] = 0;
                     else
                        geoMatrix[index(i, j, k)] = 1;
                  }
                  id++;
               }
            }
         }


         print("Perform the interpolation level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }

      print("createGeoMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   void createGeoMatrix4(std::string dataNameG, double deltax_corse)
   {
      print("createGeoMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      int size = geo_nx1*geo_nx2*geo_nx3;
      geoMatrix.resize(size, 1);

      print("Perform the solid nodes: start");
      level_interp_timer->StartTimer();

      vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
      thrFilter->SetInputData(dataSetGeo);
      thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Geometry");
      thrFilter->ThresholdBetween(1, 1);
      thrFilter->Update();
      vtkSmartPointer<vtkUnstructuredGrid> ugrid = thrFilter->GetOutput();

      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

      int numberOfCells = ugrid->GetNumberOfCells();

      double x[3];
      double xMin[3];
      double xMax[3];
      int ixMin[3];
      int ixMax[3];
      vtkIdType idc = 0;

      for (int i = 0; i < numberOfCells; i++)
      {
         vtkSmartPointer<vtkIdList> plist = vtkSmartPointer<vtkIdList>::New();
         ugrid->GetCellPoints(i, plist);
         vector <double> x1;
         vector <double> x2;
         vector <double> x3;
         for (int p = 0; p < plist->GetNumberOfIds(); p++)
         {
            ugrid->GetPoint(plist->GetId(p), x);
            x1.push_back(x[0]);
            x2.push_back(x[1]);
            x3.push_back(x[2]);
         }
         xMin[0] = *min_element(x1.begin(), x1.end());
         xMin[1] = *min_element(x2.begin(), x2.end());
         xMin[2] = *min_element(x3.begin(), x3.end());

         xMax[0] = *max_element(x1.begin(), x1.end());
         xMax[1] = *max_element(x2.begin(), x2.end());
         xMax[2] = *max_element(x3.begin(), x3.end());

         getNodeIndexes(xMin, ixMin, geo_origin, deltax);
         getNodeIndexes(xMax, ixMax, geo_origin, deltax);

         for (int k = ixMin[2]; k <= ixMax[2]; k++)
         {
            for (int j = ixMin[1]; j <= ixMax[1]; j++)
            {
               for (int i = ixMin[0]; i <= ixMax[0]; i++)
               {
                  if (i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                  {
                     geoMatrix[index(i, j, k)] = 0;
                  }
               }
            }
         }
      }
      print("Perform the solid nodes: end");
      level_interp_timer->StopTimer();
      print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");

      print("createGeoMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   void createMQMatrix4(string dataNameMQ, std::string dataNameG, double deltax_corse)
   {
      print("createMQMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetMQ(ReadDataSet(dataNameMQ.c_str()));
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      vtkSmartPointer<vtkDoubleArray> lArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
      lArrayNew->SetNumberOfComponents(1);
      lArrayNew->SetName("Level");

      int numberOfPoints = dataSetMQ->GetNumberOfPoints();

      double x[3];
      double xMin[3];
      double xMax[3];
      int ixMin[3];
      int ixMax[3];

      print("add level information: start");
      timer->StartTimer();

      for (int i = 0; i < numberOfPoints; i++)
      {
         dataSetMQ->GetPoint(i, x);
         vtkIdType id = dataSetGeo->FindPoint(x);
         if (!((int)id < 0))
         {
            lArrayNew->InsertNextTuple1(lArray->GetTuple1(id));
         }
      }

      dataSetMQ->GetPointData()->AddArray(lArrayNew);

      print("add level information: end");
      timer->StopTimer();
      print("add level information time: " + toString(timer->GetElapsedTime()) + " s");


      double range[2];
      lArray->GetRange(range);

      for (int level = (int)range[0]; level <=(int)range[1]; level++)
      {
         print("Perform the interpolation level " + toString(level) + " : start");
         level_interp_timer->StartTimer();

         print("actual memory usage befor interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

         print("Create the level grid " + toString(level) + " : start");
         level_grid_timer->StartTimer();

         double delta = deltax_corse/(double)(1<<level);

         vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
         thrFilter->SetInputData(dataSetMQ);
         thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Level");
         thrFilter->ThresholdBetween(level, level);
         thrFilter->Update();
         vtkSmartPointer<vtkUnstructuredGrid> lugrid = thrFilter->GetOutput();

         vtkSmartPointer<vtkDataArray> vxArray = lugrid->GetPointData()->GetArray("Vx");
         vtkSmartPointer<vtkDataArray> vyArray = lugrid->GetPointData()->GetArray("Vy");
         vtkSmartPointer<vtkDataArray> vzArray = lugrid->GetPointData()->GetArray("Vz");
         vtkSmartPointer<vtkDataArray> prArray = lugrid->GetPointData()->GetArray("Press");

         //vtkSmartPointer<vtkXMLUnstructuredGridWriter> uwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
         //uwriter->SetInputData(lugrid);
         //string filenameugrid = "ugrid"+toString<int>(level)+".vtu";
         //uwriter->SetFileName(filenameugrid.c_str());
         //uwriter->SetDataModeToAscii();
         //uwriter->SetCompressorTypeToNone();
         //uwriter->Update();

         if (level <(int)range[1])
         {

            int numberOfCells = lugrid->GetNumberOfCells();

            vtkIdType idc = 0;
            double val_temp[8];
            double val[8];

            for (int i = 0; i < numberOfCells; i++)
            {
               vtkIdType npts = 8;
               vtkIdType* pts;
               lugrid->GetCellPoints(i, npts, pts);
               vector <double> x1;
               vector <double> x2;
               vector <double> x3;

               for (int p = 0; p < 8; p++)
               {
                  double x[3];
                  lugrid->GetPoint(pts[p], x);
                  x1.push_back(x[0]);
                  x2.push_back(x[1]);
                  x3.push_back(x[2]);
                  val_temp[p] = vxArray->GetTuple1(pts[p]);
               }

               xMin[0] = *min_element(x1.begin(), x1.end());
               xMin[1] = *min_element(x2.begin(), x2.end());
               xMin[2] = *min_element(x3.begin(), x3.end());

               xMax[0] = *max_element(x1.begin(), x1.end());
               xMax[1] = *max_element(x2.begin(), x2.end());
               xMax[2] = *max_element(x3.begin(), x3.end());

               delta = xMax[2] - xMin[2];

               for (int p = 0; p < 8; p++)
               {
                  double x[3];
                  x[0] = x1[p];
                  x[1] = x2[p];
                  x[2] = x3[p];
                  int ix[3];
                  getNodeIndexes(x, ix, xMin, delta);
                  val[2 * (2 * ix[0] + ix[1]) + ix[2]] = val_temp[p];
               }

               getNodeIndexes(xMin, ixMin, geo_origin, deltax);
               getNodeIndexes(xMax, ixMax, geo_origin, deltax);

               double xx[3];
               double yy[3];
               double zz[3];

               xx[1] = xMin[0];
               xx[2] = xMax[0];

               yy[1] = xMin[1];
               yy[2] = xMax[1];

               zz[1] = xMin[2];
               zz[2] = xMax[2];

               for (int k = ixMin[2]; k <= ixMax[2]; k++)
               {
                  for (int j = ixMin[1]; j <= ixMax[1]; j++)
                  {
                     for (int i = ixMin[0]; i <= ixMax[0]; i++)
                     {
                        if (i >= 0 && j >= 0 && k >= 0 && i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                        {
                           double xNode[3];
                           int ix[3];
                           ix[0] = i;
                           ix[1] = j;
                           ix[2] = k;
                           getNodeCoordinates(xNode, ix, geo_origin, geo_spacing[0]);
                           xx[0] = xNode[0];
                           yy[0] = xNode[1];
                           zz[0] = xNode[2];
                           if (geoMatrix[index(i, j, k)] == 1)
                           {
                              vxMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, val);
                           }
                        }
                     }
                  }
               }
            }
         }
         else
         {
            int numberOfPoints = lugrid->GetNumberOfPoints();

            for (int i = 0; i < numberOfPoints; i++)
            {
               double x[3];
               int ix[3];
               lugrid->GetPoint(i, x);
               getNodeIndexes(x, ix, geo_origin, deltax);
               if (ix[0] >= 0 && ix[1] >= 0 && ix[2] >= 0 && ix[0] < geo_nx1 && ix[1] < geo_nx2 && ix[2] < geo_nx3 && geoMatrix[index(ix[0], ix[1], ix[2])] == 1)
                  vxMatrix[index(ix[0], ix[1], ix[2])] = vxArray->GetTuple1(i);
            }
         }

         print("Perform the interpolation level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }



      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void createMQMatrix5(string dataNameMQ, std::string dataNameG, double deltax_corse)
   {
      print("createMQMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetMQ(ReadDataSet(dataNameMQ.c_str()));
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      double xMin[3];
      double xMax[3];
      int ixMin[3];
      int ixMax[3];


      print("Perform the interpolation level: start");
      level_interp_timer->StartTimer();

      print("actual memory usage before interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

      vtkSmartPointer<vtkDataArray> vxArray = dataSetMQ->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSetMQ->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSetMQ->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Press");

      vtkSmartPointer<vtkUnstructuredGrid> lugrid = vtkUnstructuredGrid::SafeDownCast(dataSetMQ);

      int numberOfCells = lugrid->GetNumberOfCells();

      vtkIdType idc = 0;
      double val_temp[8];
      vector<double> val(8, 0);

      for (int i = 0; i < numberOfCells; i++)
      {
         vtkIdType npts = 8;
         vtkIdType* pts;
         lugrid->GetCellPoints(i, npts, pts);
         vector <double> x1;
         vector <double> x2;
         vector <double> x3;

         for (int p = 0; p < 8; p++)
         {
            double x[3];
            dataSetMQ->GetPoint(pts[p], x);
            x1.push_back(x[0]);
            x2.push_back(x[1]);
            x3.push_back(x[2]);
            val_temp[p] = vxArray->GetTuple1(pts[p]);
         }

         xMin[0] = *min_element(x1.begin(), x1.end());
         xMin[1] = *min_element(x2.begin(), x2.end());
         xMin[2] = *min_element(x3.begin(), x3.end());

         xMax[0] = *max_element(x1.begin(), x1.end());
         xMax[1] = *max_element(x2.begin(), x2.end());
         xMax[2] = *max_element(x3.begin(), x3.end());

         double delta = xMax[2] - xMin[2];

         for (int p = 0; p < 8; p++)
         {
            double x[3];
            x[0] = x1[p];
            x[1] = x2[p];
            x[2] = x3[p];
            int ix[3];
            getNodeIndexes(x, ix, xMin, delta);
            val[2 * (2 * ix[0] + ix[1]) + ix[2]] = val_temp[p];
         }

         getNodeIndexes(xMin, ixMin, geo_origin, deltax);
         getNodeIndexes(xMax, ixMax, geo_origin, deltax);

         double xx[3];
         double yy[3];
         double zz[3];

         xx[1] = xMin[0];
         xx[2] = xMax[0];

         yy[1] = xMin[1];
         yy[2] = xMax[1];

         zz[1] = xMin[2];
         zz[2] = xMax[2];

         for (int k = ixMin[2]; k <= ixMax[2]; k++)
         {
            for (int j = ixMin[1]; j <= ixMax[1]; j++)
            {
               for (int i = ixMin[0]; i <= ixMax[0]; i++)
               {
                  if (i >= 0 && j >= 0 && k >= 0 && i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                  {
                     double xNode[3];
                     int ix[3];
                     ix[0] = i;
                     ix[1] = j;
                     ix[2] = k;
                     getNodeCoordinates(xNode, ix, geo_origin, geo_spacing[0]);
                     xx[0] = xNode[0];
                     yy[0] = xNode[1];
                     zz[0] = xNode[2];
                     if (geoMatrix[index(i, j, k)] == 1)
                     {
                        vxMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, &val[0]);
                     }
                  }
               }
            }
         }
      }

      print("Perform the interpolation level: end");
      level_interp_timer->StopTimer();
      print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");

      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////
   void createMQMatrix3(string dataNameMQ, std::string dataNameG, double deltax_corse)
   {
      print("createMQMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetMQ = vtkSmartPointer<vtkDataSet>::Take(ReadDataSet(dataNameMQ.c_str()));
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo = vtkSmartPointer<vtkDataSet>::Take(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");
      vtkSmartPointer<vtkDataArray> vxArray = dataSetMQ->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSetMQ->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSetMQ->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Press");

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      vtkSmartPointer<vtkImageData> vmatrix = vtkSmartPointer<vtkImageData>::New();
      vmatrix->SetOrigin(geo_origin);
      vmatrix->SetExtent(geo_extent);
      vmatrix->SetSpacing(geo_spacing);


      double range[2];
      lArray->GetRange(range);

      for (int level = (int)range[0]; level <=(int)range[1]; level++)
      {
         print("Perform the interpolation level " + toString(level) + " : start");
         level_interp_timer->StartTimer();

         print("actual memory usage befor interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

         print("Create the level grid " + toString(level) + " : start");
         level_grid_timer->StartTimer();

         vtkSmartPointer<vtkDoubleArray> vxArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vxArrayNew->SetNumberOfComponents(1);
         vxArrayNew->SetName("Vx");

         vtkSmartPointer<vtkDoubleArray> vyArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vyArrayNew->SetNumberOfComponents(1);
         vyArrayNew->SetName("Vy");

         vtkSmartPointer<vtkDoubleArray> vzArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vzArrayNew->SetNumberOfComponents(1);
         vzArrayNew->SetName("Vz");

         vtkSmartPointer<vtkDoubleArray> prArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         prArrayNew->SetNumberOfComponents(1);
         prArrayNew->SetName("Press");

         double bounds[6];
         vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
         thrFilter->SetInputData(dataSetGeo);
         thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Level");
         thrFilter->ThresholdBetween(level, level);
         thrFilter->Update();
         vtkSmartPointer<vtkUnstructuredGrid> lugrid = thrFilter->GetOutput();
         lugrid->GetBounds(bounds);

         vtkSmartPointer<vtkImageData> limage = vtkImageData::New();
         double delta = deltax_corse/(double)(1<<level);

         double liorigin[3];
         liorigin[0] = bbox[0]-0.5*delta;
         liorigin[1] = bbox[2]-0.5*delta;
         liorigin[2] = bbox[4]-0.5*delta;

         int nx1 = cint((bbox[1] - bbox[0]) / deltax+2);
         int nx2 = cint((bbox[3] - bbox[2]) / deltax+2);
         int nx3 = cint((bbox[5] - bbox[4]) / deltax+2);

         print("level matrix NX1 x NX2 x NX3 = " + toString(nx1) + " x " + toString(nx2) + " x " + toString(nx3));

         int liextent[6];
         liextent[0] = 0;
         liextent[2] = 0;
         liextent[4] = 0;
         liextent[1] = nx1-1;
         liextent[3] = nx2-1;
         liextent[5] = nx3-1;

         double lispacing[3];
         lispacing[0] = delta;
         lispacing[1] = delta;
         lispacing[2] = delta;

         limage->SetExtent(liextent);
         limage->SetOrigin(liorigin[0], liorigin[1], liorigin[2]);
         limage->SetSpacing(delta, delta, delta);

         //double liorigin[3];
         //liorigin[0] = geo_origin[0];
         //liorigin[1] = geo_origin[1];
         //liorigin[2] = geo_origin[2];

         //int nx1 = geo_nx1;
         //int nx2 = geo_nx2;
         //int nx3 = geo_nx3;

         //limage->SetExtent(geo_extent);
         //limage->SetOrigin(geo_origin);
         //limage->SetSpacing(delta, delta, delta);

         int numberOfPoints = lugrid->GetNumberOfPoints();

         int msize = limage->GetNumberOfPoints();
         vector<double> lvxMatrix(msize, 0);
         vector<double> lvyMatrix(msize, 0);
         vector<double> lvzMatrix(msize, 0);
         vector<double> lprMatrix(msize, 0);

         double x[3];
         int   ix[3];

         for (int i = 0; i < numberOfPoints; i++)
         {
            lugrid->GetPoint(i, x);
            vtkIdType id = dataSetMQ->FindPoint(x);
            if (!((int)id < 0))
            {
               getNodeIndexes(x, ix, liorigin, delta);
               //print("ix0="+toString(ix[0])+" ix1="+ toString(ix[1])+" ix2="+toString(ix[2]));
               if (ix[0] >= 0 && ix[1] >= 0 && ix[2] >= 0 && ix[0] < nx1 && ix[1] < nx2 && ix[2] < nx3)
               {
                  lvxMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = vxArray->GetTuple1(id);
                  lvyMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = vyArray->GetTuple1(id);
                  lvzMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = vzArray->GetTuple1(id);
                  lprMatrix[index(ix[0], ix[1], ix[2], nx1, nx2)] = prArray->GetTuple1(id);
               }
            }
         }


         for (int x3 = 0; x3 < nx3; x3++)
            for (int x2 = 0; x2 < nx2; x2++)
               for (int x1 = 0; x1 < nx1; x1++)
               {
                  vxArrayNew->InsertNextValue(lvxMatrix[index(x1, x2, x3, nx1, nx2)]);
                  vyArrayNew->InsertNextValue(lvyMatrix[index(x1, x2, x3, nx1, nx2)]);
                  vzArrayNew->InsertNextValue(lvzMatrix[index(x1, x2, x3, nx1, nx2)]);
                  prArrayNew->InsertNextValue(lprMatrix[index(x1, x2, x3, nx1, nx2)]);
               }


         limage->GetPointData()->AddArray(vxArrayNew);
         limage->GetPointData()->AddArray(vyArrayNew);
         limage->GetPointData()->AddArray(vzArrayNew);
         limage->GetPointData()->AddArray(prArrayNew);

         vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
         writer->SetInputData(limage);
         //string lstr = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/limage"+toString(level)+".vti";
         string lstr = "d:/temp/sfb1/limage"+toString(level)+".vti";
         writer->SetFileName(lstr.c_str());
         writer->SetDataModeToAscii();
         writer->Update();

         print("Create the level grid " + toString(level) + " : end");
         level_grid_timer->StopTimer();
         print("creating time: " + toString(level_grid_timer->GetElapsedTime()) + " s");

         //continue;

         // Perform the interpolation
         vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
         probeFilter->SetSourceData(limage);

         // Interpolate 'Source' at these points
#if VTK_MAJOR_VERSION <= 5
         probeFilter->SetInput(vmatrix);
#else
         probeFilter->SetInputData(vmatrix);
#endif

         print("start interpolation");
         timer->StartTimer();

         probeFilter->Update();
         timer->StopTimer();
         print("end interpolation");
         print("interpolation time: " + toString(timer->GetElapsedTime()) + " s");

         int id = 0;

         vtkSmartPointer<vtkDataArray> vxArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vx");
         vtkSmartPointer<vtkDataArray> vyArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vy");
         vtkSmartPointer<vtkDataArray> vzArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Vz");
         vtkSmartPointer<vtkDataArray> prArrayP = probeFilter->GetOutput()->GetPointData()->GetArray("Press");


         for (int k = 0; k < geo_nx3; k++)
         {
            for (int j = 0; j < geo_nx2; j++)
            {
               for (int i = 0; i < geo_nx1; i++)
               {
                  if (i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
                  {
                     double vx = vxArrayP->GetTuple1(id);
                     if (vx != 0 && geoMatrix[index(i, j, k)] == 1 && vxMatrix[index(i, j, k)] == 0 && lvlMatrix[index(i, j, k)] == level)
                     {
                        vxMatrix[index(i, j, k)] = vx;
                        vyMatrix[index(i, j, k)] = vyArrayP->GetTuple1(id);
                        vzMatrix[index(i, j, k)] = vzArrayP->GetTuple1(id);
                        prMatrix[index(i, j, k)] = prArrayP->GetTuple1(id);

                     }
                  }
                  id++;
               }
            }
         }


         print("Perform the interpolation level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }

      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void createMatrix(std::string dataName, std::string dataNameG)
   {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataName + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSet(ReadDataSet(dataName.c_str()));
      print("read data set from " + dataName + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("read data set from " + dataNameG + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      print("read data set from " + dataNameG + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");
      vtkSmartPointer<vtkDataArray> vxArray = dataSet->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSet->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSet->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSet->GetPointData()->GetArray("Press");

      double bounds[6];
      dataSet->GetBounds(bounds);

      double origin[3];
      origin[0] = bounds[0] + deltax / 2.0;
      origin[1] = bounds[2] + deltax / 2.0;
      origin[2] = bounds[4] + deltax / 2.0;

      int extent[6];

      geo_nx1 = int((bounds[1] - bounds[0]) / deltax);
      geo_nx2 = int((bounds[3] - bounds[2]) / deltax);
      geo_nx3 = int((bounds[5] - bounds[4]) / deltax);

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      int size = geo_nx1*geo_nx2*geo_nx3;
      //geoMatrix = new int[size];
      //vxMatrix = new double[size];
      //vyMatrix = new double[size];
      //vzMatrix = new double[size];
      //prMatrix = new double[size];
      //vaVxMatrix = new double[size];
      //vaVyMatrix = new double[size];
      //vaVzMatrix = new double[size];
      //vaPrMatrix = new double[size];

      vtkSmartPointer<vtkImageData> vmatrix = vtkImageData::New();
      double delta = deltax;
      extent[0] = 0;
      extent[2] = 0;
      extent[4] = 0;
      extent[1] = geo_nx1 - 1;
      extent[3] = geo_nx2 - 1;
      extent[5] = geo_nx3 - 1;
      vmatrix->SetExtent(extent);
      vmatrix->SetOrigin(origin[0], origin[1], origin[2]);
      vmatrix->SetSpacing(delta, delta, delta);

      vtkSmartPointer<vtkDoubleArray> vxArrayMatrix = vtkSmartPointer<vtkDoubleArray>::New();
      vxArrayMatrix->SetNumberOfComponents(1);
      vxArrayMatrix->SetName("Vx");
      vxArrayMatrix->SetNumberOfTuples(vmatrix->GetNumberOfPoints());

      vtkSmartPointer<vtkDoubleArray> vyArrayMatrix = vtkSmartPointer<vtkDoubleArray>::New();
      vyArrayMatrix->SetNumberOfComponents(1);
      vyArrayMatrix->SetName("Vy");
      vyArrayMatrix->SetNumberOfTuples(vmatrix->GetNumberOfPoints());

      vtkSmartPointer<vtkDoubleArray> vzArrayMatrix = vtkSmartPointer<vtkDoubleArray>::New();
      vzArrayMatrix->SetNumberOfComponents(1);
      vzArrayMatrix->SetName("Vz");
      vzArrayMatrix->SetNumberOfTuples(vmatrix->GetNumberOfPoints());

      vtkSmartPointer<vtkDoubleArray> prArrayMatrix = vtkSmartPointer<vtkDoubleArray>::New();
      prArrayMatrix->SetNumberOfComponents(1);
      prArrayMatrix->SetName("Press");
      prArrayMatrix->SetNumberOfTuples(vmatrix->GetNumberOfPoints());

      print("Perform the interpolation geo: start");
      timer->StartTimer();

      // Perform the interpolation
      vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
      probeFilter->SetSourceData(dataSetGeo);

      // Interpolate 'Source' at these points
#if VTK_MAJOR_VERSION <= 5
      probeFilter->SetInput(image);
#else
      probeFilter->SetInputData(vmatrix);
#endif

      probeFilter->Update();

      print("Perform the interpolation geo: end");
      timer->StopTimer();

      vmatrix->GetPointData()->AddArray(probeFilter->GetOutput()->GetPointData()->GetArray("Geometry"));

      double range[2];
      lArray->GetRange(range);

      printf("Level: range[0] = %d, range[1] = %d\n", (int)range[0], (int)range[1]);

      for (int level = (int)range[0]; level <= (int)range[1]; level++)
      {
         print("Perform the interpolation level " + toString(level) + " : start");
         timer->StartTimer();

         print("actual memory usage befor interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

         vtkSmartPointer<vtkDoubleArray> vxArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vxArrayNew->SetNumberOfComponents(1);
         vxArrayNew->SetName("Vx");

         vtkSmartPointer<vtkDoubleArray> vyArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vyArrayNew->SetNumberOfComponents(1);
         vyArrayNew->SetName("Vy");

         vtkSmartPointer<vtkDoubleArray> vzArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         vzArrayNew->SetNumberOfComponents(1);
         vzArrayNew->SetName("Vz");

         vtkSmartPointer<vtkDoubleArray> prArrayNew = vtkSmartPointer<vtkDoubleArray>::New();
         prArrayNew->SetNumberOfComponents(1);
         prArrayNew->SetName("Press");

         vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

         vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
         vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

         int numberOfPoints = dataSetGeo->GetNumberOfPoints();
         int numberOfCells = dataSet->GetNumberOfCells();

         double x[3];
         vtkIdType idc = 0;

         for (int i = 0; i < numberOfCells; i++)
         {
            vtkSmartPointer<vtkIdList> plist = vtkSmartPointer<vtkIdList>::New();
            dataSet->GetCellPoints(i, plist);
            dataSet->GetPoint(plist->GetId(0), x);
            vtkIdType id = dataSetGeo->FindPoint(x);
            if (!((int)id < 0))
            {
               int l = (int)lArray->GetTuple1(id);
               if (l == level)
               {
                  vtkSmartPointer<vtkIdList> plist2 = vtkSmartPointer<vtkIdList>::New();
                  for (int p = 0; p < plist->GetNumberOfIds(); p++)
                  {
                     dataSet->GetPoint(plist->GetId(p), x);
                     points->InsertNextPoint(x);
                     vxArrayNew->InsertNextTuple1(vxArray->GetTuple1(plist->GetId(p)));
                     vyArrayNew->InsertNextTuple1(vyArray->GetTuple1(plist->GetId(p)));
                     vzArrayNew->InsertNextTuple1(vzArray->GetTuple1(plist->GetId(p)));
                     prArrayNew->InsertNextTuple1(prArray->GetTuple1(plist->GetId(p)));
                     plist2->InsertNextId(idc);
                     idc++;
                  }
                  cells->InsertNextCell(plist2);

               }
            }
         }

         ugrid->SetPoints(points);
         ugrid->SetCells(VTK_VOXEL, cells);
         ugrid->GetPointData()->AddArray(vxArrayNew);
         ugrid->GetPointData()->AddArray(vyArrayNew);
         ugrid->GetPointData()->AddArray(vzArrayNew);
         ugrid->GetPointData()->AddArray(prArrayNew);

         //vtkSmartPointer<vtkXMLUnstructuredGridWriter> uwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
         //uwriter->SetInputData(ugrid);
         //string filenameugrid = "ugrid"+toString<int>(level)+".vtu";
         //uwriter->SetFileName(filenameugrid.c_str());
         //uwriter->SetDataModeToAscii();
         //uwriter->SetCompressorTypeToNone();
         //uwriter->Update();

         // Perform the interpolation
         probeFilter->SetSourceData(ugrid);
         probeFilter->Update();

         vtkSmartPointer<vtkDataSet> probeMatrix = probeFilter->GetOutput();
         vtkSmartPointer<vtkDataArray> vxArrayP = probeMatrix->GetPointData()->GetArray("Vx");
         vtkSmartPointer<vtkDataArray> vyArrayP = probeMatrix->GetPointData()->GetArray("Vy");
         vtkSmartPointer<vtkDataArray> vzArrayP = probeMatrix->GetPointData()->GetArray("Vz");
         vtkSmartPointer<vtkDataArray> prArrayP = probeMatrix->GetPointData()->GetArray("Press");
         int npoits = probeMatrix->GetNumberOfPoints();

         for (int i = 0; i < npoits; i++)
         {
            double vx = vxArrayP->GetTuple1(i);
            double vy = vyArrayP->GetTuple1(i);
            double vz = vzArrayP->GetTuple1(i);
            double pr = prArrayP->GetTuple1(i);
            if (vx != 0)
            {
               vxArrayMatrix->SetTuple1(i, vx);
               vyArrayMatrix->SetTuple1(i, vy);
               vzArrayMatrix->SetTuple1(i, vz);
               prArrayMatrix->SetTuple1(i, pr);
            }
         }

         //vxArrayNew = NULL;
         //vyArrayNew = NULL;
         //vzArrayNew = NULL;
         //prArrayNew = NULL;
         //ugrid = NULL;
         //points = NULL;
         //cells = NULL;

         print("Perform the interpolation level " + toString(level) + " : end");
         timer->StopTimer();
         print("interpolation time: " + toString(timer->GetElapsedTime()) + " s");
         print("actual memory usage after interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");
      }

      //probeFilter = NULL;

      print("actual memory usage after interpolation of all levels: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

      vmatrix->GetPointData()->AddArray(vxArrayMatrix);
      vmatrix->GetPointData()->AddArray(vyArrayMatrix);
      vmatrix->GetPointData()->AddArray(vzArrayMatrix);
      vmatrix->GetPointData()->AddArray(prArrayMatrix);

      ////string vtkfilename = "vmatrix"+toString<int>(level)+".vti";
      //string vtkfilename = "vmatrix.vti";
      //vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
      //writer->SetInputData(vmatrix);
      //writer->SetFileName(vtkfilename.c_str());
      //writer->SetDataModeToAscii();
      //writer->Update();

      vtkSmartPointer<vtkDataArray> geoArray = vmatrix->GetPointData()->GetArray("Geometry");

      print("create matrix: start");
      timer->StartTimer();

      int id = 0;

      for (int k = 0; k < geo_nx3; k++)
      {
         for (int j = 0; j < geo_nx2; j++)
         {
            for (int i = 0; i < geo_nx1; i++)
            {
               int g = (int)geoArray->GetTuple1(id);
               if (g)
                  geoMatrix[index(i, j, k)] = 0;
               else
                  geoMatrix[index(i, j, k)] = 1;
               vxMatrix[index(i, j, k)] = vxArrayMatrix->GetTuple1(id);
               vyMatrix[index(i, j, k)] = vyArrayMatrix->GetTuple1(id);
               vzMatrix[index(i, j, k)] = vzArrayMatrix->GetTuple1(id);
               prMatrix[index(i, j, k)] = prArrayMatrix->GetTuple1(id);
               id++;
            }
         }
      }

      //vxArrayMatrix = NULL;
      //vyArrayMatrix = NULL;
      //vzArrayMatrix = NULL;
      //prArrayMatrix = NULL;

      //string vtkfilename2 = "VoxelMatrix.vti";
      ////vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
      //writer->SetInputData(image);
      //writer->SetFileName(vtkfilename.c_str());
      //writer->SetDataModeToAscii();
      //writer->Update();

      print("create matrix: end");
      timer->StopTimer();
      print("create matrix time: " + toString(timer->GetElapsedTime()) + " s");
      timer_total->StopTimer();
      print("create matrix total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void writeMatrixToImageFile(std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

      std::string vtkfilename = output + ".vti";

      print("write data set to " + vtkfilename + ": start");
      timer_write->StartTimer();

      vtkSmartPointer<vtkImageData> image = vtkImageData::New();
      double delta = deltax;
      int extent[6];
      extent[0] = 0;
      extent[2] = 0;
      extent[4] = 0;
      extent[1] = geo_nx1 - 1;
      extent[3] = geo_nx2 - 1;
      extent[5] = geo_nx3 - 1;
      image->SetExtent(extent);
      image->SetOrigin(geo_origin[0], geo_origin[1], geo_origin[2]);
      image->SetSpacing(delta, delta, delta);

      vtkSmartPointer<vtkDoubleArray> vxArray = vtkDoubleArray::New();
      vxArray->SetNumberOfComponents(1);
      vxArray->SetName("vx");

      vtkSmartPointer<vtkDoubleArray> vyArray = vtkDoubleArray::New();
      vyArray->SetNumberOfComponents(1);
      vyArray->SetName("vy");

      vtkSmartPointer<vtkDoubleArray> vzArray = vtkDoubleArray::New();
      vzArray->SetNumberOfComponents(1);
      vzArray->SetName("vz");

      vtkSmartPointer<vtkDoubleArray> prArray = vtkDoubleArray::New();
      prArray->SetNumberOfComponents(1);
      prArray->SetName("pr");

      vtkSmartPointer<vtkIntArray> geoArray = vtkIntArray::New();
      geoArray->SetNumberOfComponents(1);
      geoArray->SetName("geo");

      int size = geo_nx1*geo_nx2*geo_nx3;

      vxArray->SetArray(&vaVxMatrix[0], size, 1);
      vyArray->SetArray(&vaVyMatrix[0], size, 1);
      vzArray->SetArray(&vaVzMatrix[0], size, 1);
      prArray->SetArray(&vaPrMatrix[0], size, 1);
      geoArray->SetArray(&geoMatrix[0], size, 1);

      image->GetPointData()->AddArray(vxArray);
      image->GetPointData()->AddArray(vyArray);
      image->GetPointData()->AddArray(vzArray);
      image->GetPointData()->AddArray(prArray);
      image->GetPointData()->AddArray(geoArray);

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      writer->SetDataModeToAscii();
      writer->Update();

      print("write data set: end");
      timer_write->StopTimer();
      print("write data set time: " + toString(timer_write->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void writeGeoMatrixToImageFile(std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

      std::string vtkfilename = output + ".vti";

      print("write data set to " + vtkfilename + ": start");
      timer_write->StartTimer();

      vtkSmartPointer<vtkImageData> image = vtkImageData::New();

      image->SetExtent(geo_extent);
      image->SetOrigin(geo_origin);
      image->SetSpacing(geo_spacing);

      vtkSmartPointer<vtkIntArray> geoArray = vtkIntArray::New();
      geoArray->SetNumberOfComponents(1);
      geoArray->SetName("geo");

      vtkSmartPointer<vtkIntArray> lvlArray = vtkIntArray::New();
      lvlArray->SetNumberOfComponents(1);
      lvlArray->SetName("level");

      int size = geo_nx1*geo_nx2*geo_nx3;

      geoArray->SetArray(&geoMatrix[0], size, 1);
      image->GetPointData()->AddArray(geoArray);

      lvlArray->SetArray(&lvlMatrix[0], size, 1);
      image->GetPointData()->AddArray(lvlArray);

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      writer->SetDataModeToAscii();
      writer->Update();

      print("write data set: end");
      timer_write->StopTimer();
      print("write data set time: " + toString(timer_write->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void writeMQMatrixToImageFile(std::string output)
   {
      vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

      std::string vtkfilename = output + ".vti";

      print("write data set to " + vtkfilename + ": start");
      timer_write->StartTimer();

      vtkSmartPointer<vtkImageData> image = vtkImageData::New();

      image->SetExtent(geo_extent);
      image->SetOrigin(geo_origin);
      image->SetSpacing(geo_spacing);

      vtkSmartPointer<vtkDoubleArray> vxArray = vtkDoubleArray::New();
      vxArray->SetNumberOfComponents(1);
      vxArray->SetName("Vx");

      vtkSmartPointer<vtkDoubleArray> vyArray = vtkDoubleArray::New();
      vyArray->SetNumberOfComponents(1);
      vyArray->SetName("Vy");

      vtkSmartPointer<vtkDoubleArray> vzArray = vtkDoubleArray::New();
      vzArray->SetNumberOfComponents(1);
      vzArray->SetName("Vz");

      vtkSmartPointer<vtkDoubleArray> prArray = vtkDoubleArray::New();
      prArray->SetNumberOfComponents(1);
      prArray->SetName("Press");

      int size = geo_nx1*geo_nx2*geo_nx3;

      vxArray->SetArray(&vxMatrix[0], size, 1);
      vyArray->SetArray(&vyMatrix[0], size, 1);
      vzArray->SetArray(&vzMatrix[0], size, 1);
      prArray->SetArray(&prMatrix[0], size, 1);

      image->GetPointData()->AddArray(vxArray);
      image->GetPointData()->AddArray(vyArray);
      image->GetPointData()->AddArray(vzArray);
      image->GetPointData()->AddArray(prArray);

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      writer->SetDataModeToAscii();
      writer->Update();

      print("write data set: end");
      timer_write->StopTimer();
      print("write data set time: " + toString(timer_write->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   double G(double x)
   {
      if (fabs(x) <= l)
         return l - fabs(x);
      else
         return 0.0;
   }
   //////////////////////////////////////////////////////////////////////////
   double m(double x1, double x2, double x3)
   {
      return (G(x1)*G(x2)*G(x3)) / lNorm;
   }
   //////////////////////////////////////////////////////////////////////////
   void correctIndex(int& x1, int& x2, int& x3)
   {
      if (x1 < 0)   x1 = geo_nx1 + x1;
      if (x1 >= geo_nx1) x1 = x1 - geo_nx2;

      if (x2 < 0)   x2 = geo_nx2 + x2;
      if (x2 >= geo_nx2) x2 = x2 - geo_nx2;

      if (x3 < 0)   x3 = 0;
      if (x3 >= geo_nx3) x3 = geo_nx3 - 1;
   }
   //////////////////////////////////////////////////////////////////////////
   void averaging()
   {
      vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();

      print("volume averaging: start");
      timer_averaging->StartTimer();

      int i = 0;

      timer_inloop->StartTimer();
      int p = 1000000;

#pragma omp parallel num_threads(8) //private(xx, yy, zz, vx)
      {
#pragma omp for 
         for (int x3 = 0; x3 < geo_nx3; x3++)
            for (int x2 = 0; x2 < geo_nx2; x2++)
               for (int x1 = 0; x1 < geo_nx1; x1++)
               {
                  if (i%p == 0 && i != 0)
                  {
                     timer_inloop->StartTimer();
                     print("point id = " + toString(i));
                  }

                  double vx = 0.0;
                  double vy = 0.0;
                  double vz = 0.0;
                  double pr = 0.0;

                  for (int z = (int)-l; z < (int)+l; z++)
                     for (int y = (int)-l; y < (int)+l; y++)
                        for (int x = (int)-l; x < (int)+l; x++)
                        {
                           int xx = x1 + x;
                           int yy = x2 + y;
                           int zz = x3 + z;
                           if (xx < 0 || xx >= geo_nx1 || yy < 0 || yy >= geo_nx2)
                           {
                              vaVxMatrix[index(x1, x2, x3)] += 0.0;
                           }
                           else
                           {
                              correctIndex(xx, yy, zz);
                              double mm = m((double)x, (double)y, (double)z);
                              double gamma = (double)geoMatrix[index(xx, yy, zz)];
                              vx += gamma*mm*vxMatrix[index(xx, yy, zz)];
                              vy += gamma*mm*vyMatrix[index(xx, yy, zz)];
                              vz += gamma*mm*vzMatrix[index(xx, yy, zz)];
                              pr += gamma*mm*prMatrix[index(xx, yy, zz)];
                           }
                        }

                  vaVxMatrix[index(x1, x2, x3)] = vx;
                  vaVyMatrix[index(x1, x2, x3)] = vy;
                  vaVzMatrix[index(x1, x2, x3)] = vz;
                  vaPrMatrix[index(x1, x2, x3)] = pr;

                  if (i%p == 0 && i != 0)
                  {
                     timer_inloop->StopTimer();
                     print("volume dx3 = " + toString(dx3));
                     print("time per " + toString(p) + " points: " + toString(timer_inloop->GetElapsedTime()) + " s");
                     print("actual memory usage: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");
                     timer_inloop->StartTimer();
                  }
                  i++;
               }

      }
      timer_averaging->StopTimer();
      print("volume averaging: end");
      print("volume averaging time: " + toString(timer_averaging->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void volumeAveragingWithProbe()
   {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      timer->StartTimer();

      //////////////////////////////////////////////////////////////////////////
      //bKanal setup
      //string dataName = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/steps/step_488750.pvtu";
      //string dataNameG = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/grid/nodes_0.pvtu";
      //string dataNameB = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/grid/blocks_0_0.bin.vtu";
      //string output = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/";

      //string dataNameGM = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/geomatrix.vti";

      //deltax = 0.665557;
      //double dxcorse = 2.66223;
      //double l_real = 20.0;
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //block_test setup
      string dataName = "d:/Data/bkanal/1/mq.vtu";
      string dataNameG = "d:/Data/bkanal/1/nodes.vtu";
      string dataNameB = "d:/Data/bkanal/1/blocks.vtu";
      string output = "d:/temp/sfb1/";
      string outputWindow = "d:/temp/sfb1/OutputWindow.txt";
      string dataNameGM = "d:/temp/sfb1/geomatrix.vti";

      deltax = 0.665557;
      double dxcorse = 2.66223;
      double l_real = 2.4;

      vtkSmartPointer<vtkFileOutputWindow> fileOutputWindow = vtkSmartPointer<vtkFileOutputWindow>::New();
      fileOutputWindow->SetFileName(outputWindow.c_str());
      fileOutputWindow->SetFlush(0);
      fileOutputWindow->SetInstance(fileOutputWindow);

      //////////////////////////////////////////////////////////////////////////

      //vtkSmartPointer<vtkDataSet> dataSetBlocks(ReadDataSet(dataNameB.c_str()));
      //dataSetBlocks->GetBounds(bbox);
      //print("Bounding Box: NX1 x NX2 x NX3 = " + toString(bbox[1]) + " x " + toString(bbox[3]) + " x " + toString(bbox[5]));

      initBoundingBox(dataNameB);

      createGeoMatrix4(dataNameG, dxcorse);
      writeGeoMatrixToImageFile(output + "geomatrix");
      readGeoMatrix(dataNameGM);

      createMQMatrix4(dataName, dataNameG, dxcorse);
      writeMQMatrixToImageFile(output + "mq_pr");


      //l = round(l_real / deltax);
      //print("l = " + toString(l));

      //dx3 = deltax*deltax*deltax;
      //lQuadrat = l*l;
      //lNorm = lQuadrat*lQuadrat*lQuadrat;

      //createMatrix(dataName, dataNameG);
      //print("actual memory usage after createMatrix(): " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

      //averaging();
      //print("actual memory usage after averaging(): " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

      //writeMatrixToImageFile(output + "vmatrixR");

      //timer->StopTimer();
      //print("total time: " + toString(timer->GetElapsedTime()) + " s");
   }

}
