#include "Averaging.h"
#include "UbLogger.h"

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
void Averaging::createGeoMatrix(std::string dataNameG, double deltax, double geo_origin[3])
{
   UBLOG(logINFO,"createGeoMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO,"read data set from " << dataNameG << ": start");
   timer->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
   UBLOG(logINFO,"read data set from " + dataNameG + ": end");
   timer->StopTimer();
   UBLOG(logINFO,"read data set time: " << toString(timer->GetElapsedTime()) + " s");

   geoMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 1);

   UBLOG(logINFO,"Perform the solid nodes: start");
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
               if (i < dimensions[0] && j < dimensions[1] && k < dimensions[2])
               {
                  geoMatrix(i, j, k) = 0;
               }
            }
         }
      }
   }
   UBLOG(logINFO,"Perform the solid nodes: end");
   level_interp_timer->StopTimer();
   UBLOG(logINFO,"interpolation time: " << toString(level_interp_timer->GetElapsedTime()) << " s");

   UBLOG(logINFO,"createGeoMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO,"total time: " << toString(timer_total->GetElapsedTime()) << " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::writeGeoMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3])
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output + ".vti";

   UBLOG(logINFO, "write data set to "<< vtkfilename << ": start");
   timer_write->StartTimer();

   vtkSmartPointer<vtkImageData> image = vtkImageData::New();

   image->SetExtent(geo_extent);
   image->SetOrigin(geo_origin);
   image->SetSpacing(geo_spacing);

   vtkSmartPointer<vtkIntArray> geoArray = vtkIntArray::New();
   geoArray->SetNumberOfComponents(1);
   geoArray->SetName("geo");

   int size = dimensions[0]*dimensions[1]*dimensions[2];

   geoArray->SetArray(geoMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);
   image->GetPointData()->AddArray(geoArray);

   vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
   writer->SetInputData(image);
   writer->SetFileName(vtkfilename.c_str());
   writer->SetDataModeToAscii();
   //writer->SetDataModeToAppended();
   //writer->SetCompressorTypeToZLib();
   writer->Update();

   UBLOG(logINFO,"write data set: end");
   timer_write->StopTimer();
   UBLOG(logINFO,"write data set time: " << toString(timer_write->GetElapsedTime()) << " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::getNodeIndexes(double x[3], int ix[3], double origin[3], double deltax)
{
   //ix[0] = cint((x[0]-origin[0])/deltax);
   //ix[1] = cint((x[1]-origin[1])/deltax);
   //ix[2] = cint((x[2]-origin[2])/deltax);

   ix[0] = (int)round((x[0] - origin[0]) / deltax);
   ix[1] = (int)round((x[1] - origin[1]) / deltax);
   ix[2] = (int)round((x[2] - origin[2]) / deltax);
}
//////////////////////////////////////////////////////////////////////////
void Averaging::createMQMatrix(std::string dataNameMQ, double deltax, double geo_origin[3])
{
   UBLOG(logINFO,"createMQMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO,"read data set from " + dataNameMQ + ": start");
   timer->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSetMQ(ReadDataSet(dataNameMQ.c_str()));
   UBLOG(logINFO,"read data set from " + dataNameMQ + ": end");
   timer->StopTimer();
   UBLOG(logINFO,"read data set time: " + toString(timer->GetElapsedTime()) + " s");

   UBLOG(logINFO,"NX1 x NX2 x NX3 = " + toString(dimensions[0]) + " x " + toString(dimensions[1]) + " x " + toString(dimensions[2]));

   int size = dimensions[0]*dimensions[1]*dimensions[2];
   vxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   prMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkUnstructuredGrid::SafeDownCast(dataSetMQ);

   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
   vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

   int numberOfCells = ugrid->GetNumberOfCells();

   vtkSmartPointer<vtkDataArray> vxArray = ugrid->GetPointData()->GetArray("Vx");
   vtkSmartPointer<vtkDataArray> vyArray = ugrid->GetPointData()->GetArray("Vy");
   vtkSmartPointer<vtkDataArray> vzArray = ugrid->GetPointData()->GetArray("Vz");
   vtkSmartPointer<vtkDataArray> prArray = ugrid->GetPointData()->GetArray("Rho");

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
      int c =0;
      for (int k = ixMin[2]; k <= ixMax[2]; k++)
      {
         for (int j = ixMin[1]; j <= ixMax[1]; j++)
         {
            for (int i = ixMin[0]; i <= ixMax[0]; i++)
            {
               if (i >=0 && i < dimensions[0] && j >=0 && j < dimensions[1] && k >=0 && k < dimensions[2])
               {
                  if (geoMatrix(i, j, k) == 1)
                  {
                     vxMatrix(i, j, k) = vxArray->GetTuple1(plist->GetId(c));
                     vyMatrix(i, j, k) = vyArray->GetTuple1(plist->GetId(c));
                     vzMatrix(i, j, k) = vzArray->GetTuple1(plist->GetId(c));
                     prMatrix(i, j, k) = prArray->GetTuple1(plist->GetId(c));
                     c++;
                  }
               }
            }
         }
      }
   }

   UBLOG(logINFO,"createMQMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO,"total time: " + toString(timer_total->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::writeMQMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3])
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output + ".vti";

   UBLOG(logINFO,"write data set to " + vtkfilename + ": start");
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

   int size = dimensions[0]*dimensions[1]*dimensions[2];

   vxArray->SetArray(vxMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);
   vyArray->SetArray(vyMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);
   vzArray->SetArray(vzMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);
   prArray->SetArray(prMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);

   image->GetPointData()->AddArray(vxArray);
   image->GetPointData()->AddArray(vyArray);
   image->GetPointData()->AddArray(vzArray);
   image->GetPointData()->AddArray(prArray);

   vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
   writer->SetInputData(image);
   writer->SetFileName(vtkfilename.c_str());
   writer->SetDataModeToAscii();
   //writer->SetDataModeToAppended();
   //writer->SetCompressorTypeToZLib();
   writer->Update();

   writer->Delete();
   image->Delete();
   vxArray->Delete();
   vyArray->Delete();
   vzArray->Delete();

   UBLOG(logINFO,"write data set: end");
   timer_write->StopTimer();
   UBLOG(logINFO,"write data set time: " + toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
double Averaging::G(double x)
{
   //if (fabs(x) <= l)
   //   return l - fabs(x);
   //else
      return 0.0;
}
//////////////////////////////////////////////////////////////////////////