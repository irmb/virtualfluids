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
#include <omp.h>
#include <mpi.h>
using namespace std;

namespace VAI
{
   double deltax;
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
   vector<int> flgMatrix;
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
      //ix[0] = cint((x[0]-origin[0])/deltax);
      //ix[1] = cint((x[1]-origin[1])/deltax);
      //ix[2] = cint((x[2]-origin[2])/deltax);

      ix[0] = (int)round((x[0]-origin[0])/deltax);
      ix[1] = (int)round((x[1]-origin[1])/deltax);
      ix[2] = (int)round((x[2]-origin[2])/deltax);
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
   //double trilinearInterpolation(double x[3], double y[3], double z[3], double val[8])
   //{
   //   double c = (x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1]);

   //   return (x[2]-x[0])*(y[2]-y[0])*(z[2]-z[0])/c*val[0]+
   //      (x[2]-x[0])*(y[2]-y[0])*(z[0]-z[1])/c*val[1]+
   //      (x[2]-x[0])*(y[0]-y[1])*(z[2]-z[0])/c*val[2]+
   //      (x[2]-x[0])*(y[0]-y[1])*(z[0]-z[1])/c*val[3]+
   //      (x[0]-x[1])*(y[2]-y[0])*(z[2]-z[0])/c*val[4]+
   //      (x[0]-x[1])*(y[2]-y[0])*(z[0]-z[1])/c*val[5]+
   //      (x[0]-x[1])*(y[0]-y[1])*(z[2]-z[0])/c*val[6]+
   //      (x[0]-x[1])*(y[0]-y[1])*(z[0]-z[1])/c*val[7];
   //}
   //////////////////////////////////////////////////////////////////////////
   void initBoundingBox(string blocksName)
   {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();

      print("read data set from " + blocksName + ": start");
      timer->StartTimer();
      vtkSmartPointer<vtkDataSet> dataSetBlocks(ReadDataSet(blocksName.c_str()));
      print("read data set from " + blocksName + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      dataSetBlocks->GetBounds(bbox);
      print("Bounding Box: NX1 x NX2 x NX3 = " + toString(bbox[1]) + " x " + toString(bbox[3]) + " x " + toString(bbox[5]));

      geo_origin[0] = bbox[0]+0.5*deltax+deltax;
      geo_origin[1] = bbox[2]+0.5*deltax+deltax;
      geo_origin[2] = bbox[4]+0.5*deltax+deltax;

      geo_nx1 = cint((bbox[1] - bbox[0]) / deltax - deltax);
      geo_nx2 = cint((bbox[3] - bbox[2]) / deltax - deltax);
      geo_nx3 = cint((bbox[5] - bbox[4]) / deltax - deltax);

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

      //int size = geo_nx1*geo_nx2*geo_nx3;
      //lvlMatrix.resize(size, 0);
      ////geoMatrix.resize(size, 0);

      //vtkSmartPointer<vtkDataArray> lArray = dataSetBlocks->GetCellData()->GetArray("level");

      //double range[2];
      //lArray->GetRange(range);

      //for (int level = (int)range[0]; level <=(int)range[1]; level++)
      //{
      //   print("Perform the blocks level " + toString(level) + " : start");
      //   level_interp_timer->StartTimer();

      //   vtkSmartPointer<vtkThreshold> thrFilter = vtkSmartPointer<vtkThreshold>::New();
      //   thrFilter->SetInputData(dataSetBlocks);
      //   thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "level");
      //   thrFilter->ThresholdBetween(level, level);
      //   thrFilter->Update();
      //   vtkSmartPointer<vtkUnstructuredGrid> lugrid = thrFilter->GetOutput();

      //   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      //   vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

      //   int numberOfCells = lugrid->GetNumberOfCells();

      //   double x[3];
      //   double xMin[3];
      //   double xMax[3];
      //   int ixMin[3];
      //   int ixMax[3];
      //   vtkIdType idc = 0;

      //   for (int i = 0; i < numberOfCells; i++)
      //   {
      //      vtkIdType npts = 8;
      //      vtkIdType* pts;
      //      lugrid->GetCellPoints(i, npts, pts);
      //      vector <double> x1;
      //      vector <double> x2;
      //      vector <double> x3;
      //      for (int p = 0; p < 8; p++)
      //      {
      //         lugrid->GetPoint(pts[p], x);
      //         x1.push_back(x[0]);
      //         x2.push_back(x[1]);
      //         x3.push_back(x[2]);
      //      }
      //      xMin[0] = *min_element(x1.begin(), x1.end());
      //      xMin[1] = *min_element(x2.begin(), x2.end());
      //      xMin[2] = *min_element(x3.begin(), x3.end());

      //      xMax[0] = *max_element(x1.begin(), x1.end());
      //      xMax[1] = *max_element(x2.begin(), x2.end());
      //      xMax[2] = *max_element(x3.begin(), x3.end());

      //      getNodeIndexes(xMin, ixMin, geo_origin, deltax);
      //      getNodeIndexes(xMax, ixMax, geo_origin, deltax);

      //      for (int k = ixMin[2]; k <= ixMax[2]; k++)
      //      {
      //         for (int j = ixMin[1]; j <= ixMax[1]; j++)
      //         {
      //            for (int i = ixMin[0]; i <= ixMax[0]; i++)
      //            {
      //               if (i >= 0 && j >= 0 && k >= 0 && i < geo_nx1 && j < geo_nx2 && k < geo_nx3)
      //               {
      //                  lvlMatrix[index(i, j, k)] = level;
      //               }
      //            }
      //         }
      //      }
      //   }
      //   print("Perform the blocks level " + toString(level) + " : end");
      //   level_interp_timer->StopTimer();
      //   print("time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      //}
   }
   //////////////////////////////////////////////////////////////////////////

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
      //lvlMatrix.resize(size, 0);

      vtkSmartPointer<vtkDataArray> geoArray = dataSetGeo->GetPointData()->GetArray("geo");
      //vtkSmartPointer<vtkDataArray> lvlArray = dataSetGeo->GetPointData()->GetArray("level");

      int numberOfPoints = dataSetGeo->GetNumberOfPoints();

      for (int i = 0; i < numberOfPoints; i++)
      {
         geoMatrix[i] = (int)geoArray->GetTuple1(i);
         //lvlMatrix[i] = (int)lvlArray->GetTuple1(i);
      }

      print("readGeoMatrix:end");
   }
   //////////////////////////////////////////////////////////////////////////
   void readMQMatrix(string dataNameMQ)
   {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();

      print("readMQMatrix:start");

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkDataSet* dataSetMQ = ReadDataSet(dataNameMQ.c_str());
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      vtkSmartPointer<vtkImageData> image = vtkImageData::SafeDownCast(dataSetMQ);

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      vtkSmartPointer<vtkDataArray> vxArray = dataSetMQ->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSetMQ->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSetMQ->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Press");

      int numberOfPoints = dataSetMQ->GetNumberOfPoints();

      for (int i = 0; i < numberOfPoints; i++)
      {
         vxMatrix[i] = vxArray->GetTuple1(i);
         vyMatrix[i] = vyArray->GetTuple1(i);
         vzMatrix[i] = vzArray->GetTuple1(i);
         prMatrix[i] = prArray->GetTuple1(i);
      }

      dataSetMQ->Delete();
      print("readMQMatrix:end");
   }
   //////////////////////////////////////////////////////////////////////////
   void createGeoMatrix(std::string dataNameG, double deltax_corse)
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
   void createMQMatrix(string dataNameMQ, std::string dataNameG, double deltax_corse)
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

      flgMatrix.resize(size, 0);

      double xMin[3];
      double xMax[3];
      int ixMin[3];
      int ixMax[3];


      print("Perform the interpolation: start");
      level_interp_timer->StartTimer();

      print("actual memory usage before interpolation: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

      //vtkSmartPointer<vtkDataArray> vxArray = dataSetMQ->GetPointData()->GetArray("Vx");
      vtkSmartPointer<vtkDataArray> vyArray = dataSetMQ->GetPointData()->GetArray("Vy");
      vtkSmartPointer<vtkDataArray> vzArray = dataSetMQ->GetPointData()->GetArray("Vz");
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Press");

      vtkSmartPointer<vtkUnstructuredGrid> lugrid = vtkUnstructuredGrid::SafeDownCast(dataSetMQ);

      int numberOfCells = lugrid->GetNumberOfCells();

      vtkSmartPointer<vtkDataArray> vxArray = lugrid->GetPointData()->GetArray("Vx");

      double val_temp[8];
      double val[8];
      double val_x[8];
      //vector<double> val(8, 0);

      int old_level = -2;
      int new_level = -1;

      vector <vtkIdType> lpts;

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
            val_x[2 * (2 * ix[2] + ix[1]) + ix[0]] = val_temp[p];
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

         int c = 0;


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
                        //if (ixMin[0]+1 == ixMax[0] && ixMin[1]+1 == ixMax[1] && ixMin[2]+1 == ixMax[2])
                        if (lvlMatrix[index(i, j, k)] == 2)
                        {
                           vxMatrix[index(i, j, k)] = val_x[c];
                           flgMatrix[index(i, j, k)] = 1;// lvlMatrix[index(i, j, k)];
                           ////print("ohne Interpolation: index " +toString(ix[0])+","+toString(ix[1])+","+ toString(ix[2]));
                           //new_level = lvlMatrix[index(i, j, k)];
                           //if (new_level != old_level )
                           //{
                           //   print("level ohne: "+toString(lvlMatrix[index(i, j, k)]));
                           //   old_level = new_level;
                           //   
                           //}
                           //lpts.insert(lpts.begin(), pts, pts+8);

                        }
                        else
                        {
                           if (flgMatrix[index(i, j, k)] == 0)
                              vxMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, val);
                           //new_level = lvlMatrix[index(i, j, k)];
                           //if (new_level != old_level)
                           //{
                           //   print("level: "+toString(lvlMatrix[index(i, j, k)]));
                           //   old_level = new_level;
                           //   
                           //}
                        }

                     }

                  }
                  c++;
               }
            }
         }
      }

      //for (int i = 0; i < lpts.size(); i++)
      //{
      //   double x[3];
      //   int ix[3];
      //   lugrid->GetPoint(lpts[i], x);
      //   getNodeIndexes(x, ix, geo_origin, deltax);
      //   if (ix[0] >= 0 && ix[1] >= 0 && ix[2] >= 0 && ix[0] < geo_nx1 && ix[1] < geo_nx2 && ix[2] < geo_nx3 && geoMatrix[index(ix[0], ix[1], ix[2])] == 1)
      //      vxMatrix[index(ix[0], ix[1], ix[2])] = vxArray->GetTuple1(lpts[i]);
      //}

      print("Perform the interpolation: end");
      level_interp_timer->StopTimer();
      print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");

      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void createMQMatrixWithLevel(string dataNameMQ, vtkDataSet* dataSetGeo/*std::string dataNameG*/, double deltax_corse)
   {
      print("createMQMatrix:start");

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
      vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

      timer_total->StartTimer();

      print("read data set from " + dataNameMQ + ": start");
      timer->StartTimer();
      vtkDataSet* dataSetMQ = ReadDataSet(dataNameMQ.c_str());
      print("read data set from " + dataNameMQ + ": end");
      timer->StopTimer();
      print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      //print("read data set from " + dataNameG + ": start");
      //timer->StartTimer();
      //vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      //print("read data set from " + dataNameG + ": end");
      //timer->StopTimer();
      //print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      vtkSmartPointer<vtkDataArray> lArray = dataSetGeo->GetPointData()->GetArray("Level");

      int size = geo_nx1*geo_nx2*geo_nx3;
      vxMatrix.resize(size, 0);
      vyMatrix.resize(size, 0);
      vzMatrix.resize(size, 0);
      prMatrix.resize(size, 0);

      vtkDoubleArray* lArrayNew = vtkDoubleArray::New();
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

         vtkThreshold* thrFilter = vtkThreshold::New();
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
            vector <double> vx_temp(8, 0);
            vector <double> vx(8, 0);
            vector <double> vy_temp(8, 0);
            vector <double> vy(8, 0);
            vector <double> vz_temp(8, 0);
            vector <double> vz(8, 0);
            vector <double> pr_temp(8, 0);
            vector <double> pr(8, 0);

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
                  vx_temp[p] = vxArray->GetTuple1(pts[p]);
                  vy_temp[p] = vyArray->GetTuple1(pts[p]);
                  vz_temp[p] = vzArray->GetTuple1(pts[p]);
                  pr_temp[p] = prArray->GetTuple1(pts[p]);
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
                  vx[2 * (2 * ix[0] + ix[1]) + ix[2]] = vx_temp[p];
                  vy[2 * (2 * ix[0] + ix[1]) + ix[2]] = vy_temp[p];
                  vz[2 * (2 * ix[0] + ix[1]) + ix[2]] = vz_temp[p];
                  pr[2 * (2 * ix[0] + ix[1]) + ix[2]] = pr_temp[p];
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
                              vxMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, &vx[0]);
                              vyMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, &vy[0]);
                              vzMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, &vz[0]);
                              prMatrix[index(i, j, k)] = trilinearInterpolation(xx, yy, zz, &pr[0]);
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
               {
                  vxMatrix[index(ix[0], ix[1], ix[2])] = vxArray->GetTuple1(i);
                  vyMatrix[index(ix[0], ix[1], ix[2])] = vyArray->GetTuple1(i);
                  vzMatrix[index(ix[0], ix[1], ix[2])] = vzArray->GetTuple1(i);
                  prMatrix[index(ix[0], ix[1], ix[2])] = prArray->GetTuple1(i);
               }
            }
         }

         thrFilter->Delete();
         print("Perform the interpolation level " + toString(level) + " : end");
         level_interp_timer->StopTimer();
         print("interpolation time: " + toString(level_interp_timer->GetElapsedTime()) + " s");
      }

      dataSetMQ->Delete();
      lArrayNew->Delete();

      print("createMQMatrix:end");
      timer_total->StopTimer();
      print("total time: " + toString(timer_total->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void writeVAMatrixToImageFile(std::string output)
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

      int size = geo_nx1*geo_nx2*geo_nx3;

      vxArray->SetArray(&vaVxMatrix[0], size, 1);
      vyArray->SetArray(&vaVyMatrix[0], size, 1);
      vzArray->SetArray(&vaVzMatrix[0], size, 1);
      prArray->SetArray(&vaPrMatrix[0], size, 1);

      image->GetPointData()->AddArray(vxArray);
      image->GetPointData()->AddArray(vyArray);
      image->GetPointData()->AddArray(vzArray);
      image->GetPointData()->AddArray(prArray);

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      //writer->SetDataModeToAscii();
      writer->SetDataModeToAppended();
      writer->SetCompressorTypeToZLib();
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

      //lvlArray->SetArray(&lvlMatrix[0], size, 1);
      //image->GetPointData()->AddArray(lvlArray);

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      //writer->SetDataModeToAscii();
      writer->SetDataModeToAppended();
      writer->SetCompressorTypeToZLib();
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

      vtkImageData* image = vtkImageData::New();

      image->SetExtent(geo_extent);
      image->SetOrigin(geo_origin);
      image->SetSpacing(geo_spacing);

      vtkDoubleArray* vxArray = vtkDoubleArray::New();
      vxArray->SetNumberOfComponents(1);
      vxArray->SetName("Vx");

      vtkDoubleArray* vyArray = vtkDoubleArray::New();
      vyArray->SetNumberOfComponents(1);
      vyArray->SetName("Vy");

      vtkDoubleArray* vzArray = vtkDoubleArray::New();
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

      vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      //writer->SetDataModeToAscii();
      writer->SetDataModeToAppended();
      //writer->SetCompressorTypeToZLib();
      writer->Update();

      writer->Delete();
      image->Delete();
      vxArray->Delete();
      vyArray->Delete();
      vzArray->Delete();

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
      if (x1 >= geo_nx1) x1 = x1 - geo_nx1;

      if (x2 < 0)   x2 = geo_nx2 + x2;
      if (x2 >= geo_nx2) x2 = x2 - geo_nx2;

      if (x3 < 0)   x3 = 0;
      if (x3 >= geo_nx3) x3 = geo_nx3 - 1;
   }
   //////////////////////////////////////////////////////////////////////////
   void averaging(double l_real)
   {
      vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

      print("volume averaging: start");
      //timer_averaging->StartTimer();

      l = round(l_real / deltax);
      print("l = " + toString(l));

      lQuadrat = l*l;
      lNorm = lQuadrat*lQuadrat*lQuadrat;

      int size = geo_nx1*geo_nx2*geo_nx3;
      vaVxMatrix.resize(size, 0);
      vaVyMatrix.resize(size, 0);
      vaVzMatrix.resize(size, 0);
      vaPrMatrix.resize(size, 0);

      vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();
      //timer_inloop->StartTimer();
      int p = 1000000;

      omp_set_num_threads(8);

      //#pragma omp parallel num_threads(4) //private(i)
      {
         int i = 0;
#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
         for (int x3 = 0; x3 < geo_nx3; x3++)
            for (int x2 = 0; x2 < geo_nx2; x2++)
               for (int x1 = 0; x1 < geo_nx1; x1++)
               {
                  //int ID = omp_get_thread_num();
                  //if (i == 0 && ID == 0)
                  //{
                  //   timer_inloop->StartTimer();
                  //   print("point id = " + toString(i));
                  //}
                  double vx = 0.0;
                  double vy = 0.0;
                  double vz = 0.0;
                  double pr = 0.0;
                  
                  int ll = (int)l;

                  //#pragma omp parallel for
                  for (int z = -ll; z <= +ll; z++)
                     for (int y = -ll; y <= +ll; y++)
                        for (int x = -ll; x <= +ll; x++)
                        {
                           int xx = x1 + x;
                           int yy = x2 + y;
                           int zz = x3 + z;

                           //correctIndex(xx, yy, zz);
                           if (xx < 0)   xx = geo_nx1 + xx;
                           if (xx >= geo_nx1) xx = xx - geo_nx1;

                           if (yy < 0)   yy = geo_nx2 + yy;
                           if (yy >= geo_nx2) yy = yy - geo_nx2;

                           if (zz < 0)   zz = 0;
                           if (zz >= geo_nx3) zz = geo_nx3 - 1;

                           int indx = index(xx, yy, zz);

                           //double mm = m((double)x, (double)y, (double)z);
                           double mm = (G((double)x)*G((double)y)*G((double)z)) / lNorm;
                           double gamma = (double)geoMatrix[indx];

                           vx += gamma*mm*vxMatrix[indx];
                           vy += gamma*mm*vyMatrix[indx];
                           vz += gamma*mm*vzMatrix[indx];
                           pr += gamma*mm*prMatrix[indx];

                        }

                  int indx = index(x1, x2, x3);
                  vaVxMatrix[indx] = vx;
                  vaVyMatrix[indx] = vy;
                  vaVzMatrix[indx] = vz;
                  vaPrMatrix[indx] = pr;

                  //if (i%p == 0 && i != 0 && ID == 0)
                  //{
                  //   timer_inloop->StopTimer();
                  //   print("point id = " + toString(i));
                  //   print("time per " + toString(p) + " points: " + toString(timer_inloop->GetElapsedTime()) + " s");
                  //   print("actual memory usage: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");
                  //   timer_inloop->StartTimer();
                  //   print("thread id: "+toString(ID));
                  //   print("Number of treads: "+toString(omp_get_num_threads()));
                  //}
                  //i++;
               }

      }
      timer_averaging->StopTimer();
      print("volume averaging: end");
      print("volume averaging time: " + toString(timer_averaging->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void averagingWithMPI(double l_real)
   {
      vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

      print("volume averaging: start");
      //timer_averaging->StartTimer();

      l = round(l_real / deltax);
      print("l = " + toString(l));

      lQuadrat = l*l;
      lNorm = lQuadrat*lQuadrat*lQuadrat;

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      int size = geo_nx1*geo_nx2*geo_nx3;
      vaVxMatrix.resize(size, 0);
      //vaVyMatrix.resize(size, 0);
      //vaVzMatrix.resize(size, 0);
      //vaPrMatrix.resize(size, 0);

      int numprocs, PID;
      MPI_Comm_rank(MPI_COMM_WORLD, &PID);
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

      int part = (int)round((double)geo_nx1 / (double)numprocs);
      print("part = " + toString(part));

      int startX1 = part * PID;
      int stopX1 = startX1 + part;
      if (PID == numprocs-1)
      {
         stopX1 = geo_nx1;
      }

      print("startX1 = " + toString(startX1));
      print("stopX1 = " + toString(stopX1));

      vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();
      //timer_inloop->StartTimer();
      int p = 1000000;

      omp_set_num_threads(8);

      //#pragma omp parallel num_threads(4) //private(i)
      {
         int i = 0;
#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
         for (int x3 = 0; x3 < geo_nx3; x3++)
            for (int x2 = 0; x2 < geo_nx2; x2++)
               for (int x1 = startX1; x1 < stopX1; x1++)
               {
                  //int ID = omp_get_thread_num();
                  //if (i == 0 && ID == 0)
                  //{
                  //   timer_inloop->StartTimer();
                  //   print("point id = " + toString(i));
                  //}
                  double vx = 0.0;
                  double vy = 0.0;
                  double vz = 0.0;
                  double pr = 0.0;

                  int ll = (int)l;

                  //#pragma omp parallel for
                  for (int z = -ll; z <= +ll; z++)
                     for (int y = -ll; y <= +ll; y++)
                        for (int x = -ll; x <= +ll; x++)
                        {
                           int xx = x1 + x;
                           int yy = x2 + y;
                           int zz = x3 + z;

                           //correctIndex(xx, yy, zz);
                           if (xx < 0)   xx = geo_nx1 + xx;
                           if (xx >= geo_nx1) xx = xx - geo_nx1;

                           if (yy < 0)   yy = geo_nx2 + yy;
                           if (yy >= geo_nx2) yy = yy - geo_nx2;

                           if (zz < 0)   zz = 0;
                           if (zz >= geo_nx3) zz = geo_nx3 - 1;

                           int indx = index(xx, yy, zz);

                           //double mm = m((double)x, (double)y, (double)z);
                           double mm = (G((double)x)*G((double)y)*G((double)z)) / lNorm;
                           double gamma = (double)geoMatrix[indx];

                           vx += gamma*mm*vxMatrix[indx];
                           //vy += gamma*mm*vyMatrix[indx];
                           //vz += gamma*mm*vzMatrix[indx];
                           //pr += gamma*mm*prMatrix[indx];

                        }

                  int indx = index(x1, x2, x3);
                  vaVxMatrix[indx] = vx;
                  //vaVyMatrix[indx] = vy;
                  //vaVzMatrix[indx] = vz;
                  //vaPrMatrix[indx] = pr;

                  //if (i%p == 0 && i != 0 && ID == 0)
                  //{
                  //   timer_inloop->StopTimer();
                  //   print("point id = " + toString(i));
                  //   print("time per " + toString(p) + " points: " + toString(timer_inloop->GetElapsedTime()) + " s");
                  //   print("actual memory usage: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");
                  //   timer_inloop->StartTimer();
                  //   print("thread id: "+toString(ID));
                  //   print("Number of treads: "+toString(omp_get_num_threads()));
                  //}
                  //i++;
               }

      }


      if (PID == 0)
      {
         vector<double> receiveBuffer;
         for (int i = 1; i < numprocs; i++)
         {
            int count, lstartX1, lstopX1;
            MPI_Status status;
            MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            receiveBuffer.resize(count);
            MPI_Recv(&receiveBuffer[0], count, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&lstartX1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&lstopX1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            int c = 0;
            for (int x3 = 0; x3 < geo_nx3; x3++)
               for (int x2 = 0; x2 < geo_nx2; x2++)
                  for (int x1 = lstartX1; x1 < lstopX1; x1++)
                  {
                     int indx = index(x1, x2, x3);
                     vaVxMatrix[indx] = receiveBuffer[c];
                     c++;
                  }
         }
      } 
      else
      {
         vector<double> sendBuffer;
         for (int x3 = 0; x3 < geo_nx3; x3++)
            for (int x2 = 0; x2 < geo_nx2; x2++)
               for (int x1 = startX1; x1 < stopX1; x1++)
               {
                  int indx = index(x1, x2, x3);
                  sendBuffer.push_back(vaVxMatrix[indx]);
               }
         int count = sendBuffer.size();
         MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
         MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
         MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
         MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }


      timer_averaging->StopTimer();
      print("volume averaging: end");
      print("volume averaging time: " + toString(timer_averaging->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
   void volumeAveragingWithInterpolation(int n)
   {
      //int required = MPI_THREAD_SERIALIZED;  // Required level of MPI threading support

      ///* Each thread will call MPI routines, but these calls will be coordinated
      //to occur only one at a time within a process. */

      //int provided; // Provided level of MPI threading support

      //MPI_Init_thread(NULL, NULL, required, &provided);

      //int numprocs, PID;
      //MPI_Comm_rank(MPI_COMM_WORLD, &PID);
      //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      //timer->StartTimer();

      //////////////////////////////////////////////////////////////////////////
      //bKanal setup
      string dataName = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/steps/step_488750.pvtu";
      string dataNameG = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/Geo/nodes.vtu";
      string dataNameB = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/grid/blocks_0_0.bin.vtu";

      //string dataName = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/test/mq.vtu";
      //string dataNameG = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/test/nodes.vtu";
      //string dataNameB = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/test/blocks.vtu";
      string output = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/";

      string outputi = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/Interpolation/";

      string dataNameGM = "/hpc3lustre/work/koskuche/SFB880/BKanalAvData/Interpolation/geomatrix.vti";

      deltax = 0.665557;
      double dxcorse = 2.66223;
      double l_real = 10.0 + 10.0; // l = df + dp
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //block_test setup
      //string dataName = "d:/Data/bkanal/1/mq.vtu";
      //string dataNameG = "d:/Data/bkanal/1/nodes.vtu";
      //string dataNameB = "d:/Data/bkanal/1/blocks.vtu";
      //string output = "d:/temp/sfb1/";
      //string outputi = "d:/temp/sfb1/";
      //string outputWindow = "d:/temp/sfb1/OutputWindow.txt";
      //string dataNameGM = "d:/temp/sfb1/geomatrix.vti";

      //deltax = 0.665557;
      //double dxcorse = 2.66223;
      //double l_real = 10.0 + 10.0; // l = df + dp

      //vtkSmartPointer<vtkFileOutputWindow> fileOutputWindow = vtkSmartPointer<vtkFileOutputWindow>::New();
      //fileOutputWindow->SetFileName(outputWindow.c_str());
      //fileOutputWindow->SetFlush(0);
      //fileOutputWindow->SetInstance(fileOutputWindow);

      //////////////////////////////////////////////////////////////////////////
      //omp_set_num_threads(8);
      //print("Number of threads: "+toString(omp_get_num_threads()));

      //initBoundingBox(dataNameB);
      //createGeoMatrix(dataNameG, dxcorse);
      //writeGeoMatrixToImageFile(outputi + "geomatrix");
      
      //readGeoMatrix(outputi + "geomatrix.vti");
      //writeMatrixToBinaryFiles < int >(geoMatrix, outputi + "geomatrix.bin");

      //print("read data set from " + dataNameG + ": start");
      //timer->StartTimer();
      //vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
      //print("read data set from " + dataNameG + ": end");
      //timer->StopTimer();
      //print("read data set time: " + toString(timer->GetElapsedTime()) + " s");

      //string path = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/steps/step_";
      //string ext = ".pvtu";

      int min, max;

      switch (n)
      {
      case 0: 
         min = 488750;
         max = 509000;
         break;
      case 1:
         min = 509000;
         max = 529250;
         break;
      case 2:
         min = 529250;
         max = 549500;
         break;
      case 3:
         min = 549500;
         max = 569750;
         break;
      case 4:
         min = 569750;
         max = 592250;
         break;
      default:
         break;
      }

      geo_nx1 = 899;
      geo_nx2 = 599;
      geo_nx3 = 599;

      //readMatrixFromBinaryFiles < int >(output+"Geo/geomatrix.raw", geoMatrix);

      //for (int i = min; i < max; i += 2250)
      //{
      //   timer->StartTimer();
      //   print("Create MQ Matrix: "+toString(i));
      //   //createMQMatrixWithLevel(path+toString(i)+ext, dataSetGeo, dxcorse);
      //   //writeMQMatrixToImageFile(outputi + "mq"+toString(i));
      //   readMQMatrix(outputi+"mq"+toString(i)+".vti");
      //   writeMatrixToBinaryFiles<double>(vxMatrix, output+"vxMatrix/vxMatrix"+toString(i)+".raw");
      //   writeMatrixToBinaryFiles<double>(vyMatrix, output+"vyMatrix/vyMatrix"+toString(i)+".raw");
      //   writeMatrixToBinaryFiles<double>(vzMatrix, output+"vzMatrix/vzMatrix"+toString(i)+".raw");
      //   timer->StopTimer();
      //   print("total time: " + toString(timer->GetElapsedTime()) + " s");
      //}


      //createMQMatrixWithLevel(dataName, dataSetGeo, dxcorse);
      //writeMQMatrixToImageFile(output + "mq");

      //readMQMatrix(outputi+"mq.vti");
      //readMQMatrix(outputi+"mq491000.vti");

      //writeMQMatrixToBinaryFiles(vxMatrix, output+"vxMatrix.bin");
      //readMatrixFromBinaryFiles(output+"matrixVx.bin", vxMatrix);

      //writeMQMatrixToImageFile(output + "mq_new");

      //averagingWithMPI(l_real);
      //////print("actual memory usage after averaging: " + toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");

      //if (PID == 0)
      //{
      //   //writeVAMatrixToImageFile(output + "vmatrixALL");
      //   writeMQMatrixToBinaryFiles(vxMatrix, output+"vaVxMatrix.bin");
      //}
      //writeVAMatrixToImageFile(output + "vmatrix_"+toString(PID));

      //timer->StopTimer();
      //print("total time: " + toString(timer->GetElapsedTime()) + " s");

      //vtkDataSet* dataSetMQ = ReadDataSet("/hpc3lustre/work/koskuche/SFB880/porplate/mq/mq14000.pvtu");
      //dataSetMQ->Delete();

      //MPI_Finalize();
   }

}
