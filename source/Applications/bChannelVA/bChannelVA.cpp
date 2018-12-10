#include <iostream>
#include <string>
#include "VirtualFluids.h"

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

using namespace std;
namespace VolAv
{
   double deltax;
   double l;
   double lQuadrat;
   double lNorm;

   int geo_nx1, geo_nx2, geo_nx3;
   //int geo_extent[6];
   //double geo_origin[3];
   //double geo_spacing[3];

   static vector<int> geoMatrix;
   vector<double> vxMatrix;
   vector<double> vyMatrix;
   vector<double> vzMatrix;
   vector<double> prMatrix;

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

      ix[0] = (int)round((x[0] - origin[0]) / deltax);
      ix[1] = (int)round((x[1] - origin[1]) / deltax);
      ix[2] = (int)round((x[2] - origin[2]) / deltax);
   }
   //////////////////////////////////////////////////////////////////////////
   void getNodeCoordinates(double x[3], int ix[3], double origin[3], double deltax)
   {
      x[0] = origin[0] + ix[0] * deltax;
      x[1] = origin[1] + ix[1] * deltax;
      x[2] = origin[2] + ix[2] * deltax;
   }
   //////////////////////////////////////////////////////////////////////////
   void createGeoMatrix(std::string dataNameG, double deltax, double geo_origin[3])
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
   /////////////////////////////////////////////////
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
      int geo_extent[6];
      double geo_origin[3];
      double geo_spacing[3];
      image->GetExtent(geo_extent);
      image->GetOrigin(geo_origin);
      image->GetSpacing(geo_spacing);

      geo_nx1 = geo_extent[1] + 1;
      geo_nx2 = geo_extent[3] + 1;
      geo_nx3 = geo_extent[5] + 1;

      print("NX1 x NX2 x NX3 = " + toString(geo_nx1) + " x " + toString(geo_nx2) + " x " + toString(geo_nx3));

      int size = geo_nx1*geo_nx2*geo_nx3;
      geoMatrix.resize(size, 0);

      vtkSmartPointer<vtkDataArray> geoArray = dataSetGeo->GetPointData()->GetArray("geo");

      int numberOfPoints = dataSetGeo->GetNumberOfPoints();

      for (int i = 0; i < numberOfPoints; i++)
      {
         geoMatrix[i] = (int)geoArray->GetTuple1(i);
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
      vtkSmartPointer<vtkDataArray> prArray = dataSetMQ->GetPointData()->GetArray("Rho");

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
   void writeGeoMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3])
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

      int size = geo_nx1*geo_nx2*geo_nx3;

      //geoArray->SetArray(&geoMatrix[0], size, 1);
      vector<int> testmat(size, 1);
      geoArray->SetArray(&testmat[0], size, 1);
      image->GetPointData()->AddArray(geoArray);

      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
      writer->SetInputData(image);
      writer->SetFileName(vtkfilename.c_str());
      writer->SetDataModeToAscii();
      //writer->SetDataModeToAppended();
      //writer->SetCompressorTypeToZLib();
      writer->Update();

      print("write data set: end");
      timer_write->StopTimer();
      print("write data set time: " + toString(timer_write->GetElapsedTime()) + " s");
   }
   //////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   VolAv::geo_nx1 = 60;
   VolAv::geo_nx2 = 40;
   VolAv::geo_nx3 = 40;
   int geo_extent[6] = { 0,60,0,40,0,40 };
   double geo_origin[3] = { 0,0,0 };
   double geo_spacing[3] = { 1,1,1 };

   VolAv::createGeoMatrix("e:/temp/BreugemChannelAnisotrop/bc/bc0.pvtu", 5, geo_origin);

   //VolAv::readMQMatrix("e:/temp/BreugemChannelAnisotrop/mq/mq100.pvtu");

   VolAv::writeGeoMatrixToImageFile("e:/temp/BreugemChannelAnisotrop/va/geoMatrix", geo_extent, geo_origin, geo_spacing);

   return 0;
}
