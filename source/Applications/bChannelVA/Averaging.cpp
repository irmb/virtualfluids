#include "Averaging.h"
#include "UbLogger.h"

#include "ReadDataSet.h"

//#include "Postprocessing.h"
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

#include <omp.h>
#include <mpi.h>

using namespace std;
void Averaging::createGeoMatrix(std::string dataNameG, double deltax, double geo_origin[3])
{
   UBLOG(logINFO, "createGeoMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO, "read data set from " << dataNameG << ": start");
   timer->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
   UBLOG(logINFO, "read data set from " + dataNameG + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " << UbSystem::toString(timer->GetElapsedTime()) + " s");

   geoMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 1);

   UBLOG(logINFO, "Perform the solid nodes: start");
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
   UBLOG(logINFO, "Perform the solid nodes: end");
   level_interp_timer->StopTimer();
   UBLOG(logINFO, "interpolation time: " << UbSystem::toString(level_interp_timer->GetElapsedTime()) << " s");

   UBLOG(logINFO, "createGeoMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO, "total time: " << UbSystem::toString(timer_total->GetElapsedTime()) << " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::writeGeoMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3])
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output + ".vti";

   UBLOG(logINFO, "write data set to " << vtkfilename << ": start");
   timer_write->StartTimer();

   vtkSmartPointer<vtkImageData> image = vtkImageData::New();

   image->SetExtent(geo_extent);
   image->SetOrigin(geo_origin);
   image->SetSpacing(geo_spacing);

   vtkSmartPointer<vtkIntArray> geoArray = vtkIntArray::New();
   geoArray->SetNumberOfComponents(1);
   geoArray->SetName("geo");

   int size = dimensions[0] * dimensions[1] * dimensions[2];

   geoArray->SetArray(geoMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);
   image->GetPointData()->AddArray(geoArray);

   vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
   writer->SetInputData(image);
   writer->SetFileName(vtkfilename.c_str());
   //writer->SetDataModeToAscii();
   writer->SetDataModeToAppended();
   writer->SetCompressorTypeToZLib();
   writer->Update();

   UBLOG(logINFO, "write data set: end");
   timer_write->StopTimer();
   UBLOG(logINFO, "write data set time: " << UbSystem::toString(timer_write->GetElapsedTime()) << " s");
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
   UBLOG(logINFO, "createMQMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO, "read data set from " + dataNameMQ + ": start");
   timer->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSetMQ(ReadDataSet(dataNameMQ.c_str()));
   UBLOG(logINFO, "read data set from " + dataNameMQ + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " + UbSystem::toString(timer->GetElapsedTime()) + " s");

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) + " x " + UbSystem::toString(dimensions[1]) + " x " + UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
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
      int c = 0;
      for (int k = ixMin[2]; k <= ixMax[2]; k++)
      {
         for (int j = ixMin[1]; j <= ixMax[1]; j++)
         {
            for (int i = ixMin[0]; i <= ixMax[0]; i++)
            {
               if (i >= 0 && i < dimensions[0] && j >= 0 && j < dimensions[1] && k >= 0 && k < dimensions[2])
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

   UBLOG(logINFO, "createMQMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO, "total time: " + UbSystem::toString(timer_total->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::writeMQMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3])
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output + ".vti";

   UBLOG(logINFO, "write data set to " + vtkfilename + ": start");
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

   int size = dimensions[0] * dimensions[1] * dimensions[2];

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
   //writer->SetDataModeToAscii();
   writer->SetDataModeToAppended();
   writer->SetCompressorTypeToZLib();
   writer->Update();

   writer->Delete();
   image->Delete();
   vxArray->Delete();
   vyArray->Delete();
   vzArray->Delete();

   UBLOG(logINFO, "write data set: end");
   timer_write->StopTimer();
   UBLOG(logINFO, "write data set time: " + UbSystem::toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::volumeAveragingWithMPI(double l_real, double deltax)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "volume averaging: start");
   //timer_averaging->StartTimer();

   double l = round(l_real / deltax);
   UBLOG(logINFO, "l = " + UbSystem::toString(l));

   double lQuadrat = l*l;
   double lNorm = lQuadrat*lQuadrat*lQuadrat;

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) << " x " + UbSystem::toString(dimensions[1]) << " x " << UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   vaVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   int numprocs, PID;
   MPI_Comm_rank(MPI_COMM_WORLD, &PID);
   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

   int part = (int)round((double)dimensions[0] / (double)numprocs);
   UBLOG(logINFO, "part = " + UbSystem::toString(part));

   int startX1 = part * PID;
   int stopX1 = startX1 + part;
   if (PID == numprocs - 1)
   {
      stopX1 = dimensions[0];
   }

   UBLOG(logINFO, "startX1 = " + UbSystem::toString(startX1));
   UBLOG(logINFO, "stopX1 = " + UbSystem::toString(stopX1));

   vtkSmartPointer<vtkTimerLog> timer_inloop = vtkSmartPointer<vtkTimerLog>::New();
   //timer_inloop->StartTimer();
   int p = 1000000;

   //omp_set_num_threads(8);

   //#pragma omp parallel num_threads(4) //private(i)
   {
      int i = 0;
//#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               //int ID = omp_get_thread_num();
               //if (i == 0 && ID == 0)
               //{
               //   timer_inloop->StartTimer();
               //   print("point id = " + UbSystem::toString(i));
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
                        if (xx < 0)   xx = dimensions[0] + xx;
                        if (xx >= dimensions[0]) xx = xx - dimensions[0];

                        if (yy < 0)   yy = dimensions[1] + yy;
                        if (yy >= dimensions[1]) yy = yy - dimensions[1];

                        if (zz < 0)   zz = 0;
                        if (zz >= dimensions[2]) zz = dimensions[2] - 1;

                        double mm = (G((double)x, l)*G((double)y, l)*G((double)z, l)) / lNorm;
                        double gamma = (double)geoMatrix(xx, yy, zz);

                        vx += gamma*mm*vxMatrix(xx, yy, zz);
                        vy += gamma*mm*vyMatrix(xx, yy, zz);
                        vz += gamma*mm*vzMatrix(xx, yy, zz);
                        pr += gamma*mm*prMatrix(xx, yy, zz);

                     }

               vaVxMatrix(x1, x2, x3) = vx;
               vaVyMatrix(x1, x2, x3) = vy;
               vaVzMatrix(x1, x2, x3) = vz;
               vaPrMatrix(x1, x2, x3) = pr;

               //if (i%p == 0 && i != 0 && ID == 0)
               //{
               //   timer_inloop->StopTimer();
               //   print("point id = " + UbSystem::toString(i));
               //   print("time per " + UbSystem::toString(p) + " points: " + UbSystem::toString(timer_inloop->GetElapsedTime()) + " s");
               //   print("actual memory usage: " + UbSystem::toString(MemoryUtil::getPhysMemUsedByMe() / 1e9) + " GByte");
               //   timer_inloop->StartTimer();
               //   print("thread id: "+UbSystem::toString(ID));
               //   print("Number of treads: "+UbSystem::toString(omp_get_num_threads()));
               //}
               //i++;
            }

   }


   //if (PID == 0)
   //{
   //   vector<double> receiveBuffer;
   //   for (int i = 1; i < numprocs; i++)
   //   {
   //      int count, lstartX1, lstopX1;
   //      MPI_Status status;
   //      MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
   //      receiveBuffer.resize(count);
   //      MPI_Recv(&receiveBuffer[0], count, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
   //      MPI_Recv(&lstartX1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
   //      MPI_Recv(&lstopX1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
   //      int c = 0;
   //      for (int x3 = 0; x3 < dimensions[2]; x3++)
   //         for (int x2 = 0; x2 < dimensions[1]; x2++)
   //            for (int x1 = lstartX1; x1 < lstopX1; x1++)
   //            {
   //               vaVxMatrix(x1, x2, x3) = receiveBuffer[c];
   //               //vaVxMatrix(x1, x2, x3) = receiveBuffer[c];
   //               //vaVxMatrix(x1, x2, x3) = receiveBuffer[c];
   //               //vaVxMatrix(x1, x2, x3) = receiveBuffer[c];
   //               c++;
   //            }
   //   }
   //}
   //else
   //{
   //   vector<double> sendBuffer;
   //   for (int x3 = 0; x3 < dimensions[2]; x3++)
   //      for (int x2 = 0; x2 < dimensions[1]; x2++)
   //         for (int x1 = startX1; x1 < stopX1; x1++)
   //         {
   //            sendBuffer.push_back(vaVxMatrix(x1, x2, x3));
   //         }
   //   int count = (int)sendBuffer.size();
   //   MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   //   MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   //   MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   //   MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   //}

   timer_averaging->StopTimer();
   UBLOG(logINFO, "volume averaging: end");
   UBLOG(logINFO, "volume averaging time: " + UbSystem::toString(timer_averaging->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::readGeoMatrix(string dataNameG)
{
   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "readGeoMatrix:start");

   UBLOG(logINFO, "read data set from " + dataNameG + ": start");
   timer->StartTimer();
   vtkSmartPointer<vtkDataSet> dataSetGeo(ReadDataSet(dataNameG.c_str()));
   UBLOG(logINFO, "read data set from " + dataNameG + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " + UbSystem::toString(timer->GetElapsedTime()) + " s");

   vtkSmartPointer<vtkImageData> image = vtkImageData::SafeDownCast(dataSetGeo);

   int geo_extent[6];
   double geo_origin[3];
   double geo_spacing[3];

   image->GetExtent(geo_extent);
   image->GetOrigin(geo_origin);
   image->GetSpacing(geo_spacing);

   int geo_nx1 = geo_extent[1] + 1;
   int geo_nx2 = geo_extent[3] + 1;
   int geo_nx3 = geo_extent[5] + 1;

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(geo_nx1) + " x " + UbSystem::toString(geo_nx2) + " x " + UbSystem::toString(geo_nx3));

   geoMatrix.resize(geo_nx1, geo_nx2, geo_nx3, 0);

   vtkSmartPointer<vtkDataArray> geoArray = dataSetGeo->GetPointData()->GetArray("geo");

   int numberOfPoints = dataSetGeo->GetNumberOfPoints();
   int* gm = geoMatrix.getStartAdressOfSortedArray(0, 0, 0);
   for (int i = 0; i < numberOfPoints; i++)
   {
      gm[i] = (int)geoArray->GetTuple1(i);
   }

   UBLOG(logINFO, "readGeoMatrix:end");
}
void Averaging::writeGeoMatrixToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<int>(geoMatrix.getDataVector(), fname);
}
void Averaging::readGeoMatrixFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<int>(fname, geoMatrix.getDataVector());
}
void Averaging::initVolumeAveragingValues()
{
   sumVaVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanVolumeAveragingValues()
{
   meanVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initFluctuationsofVolumeAveragingValues()
{
   FlucVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   FlucVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   FlucVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   FlucPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initSumOfFluctuations()
{
   sumFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanOfFluctuations()
{
   meanFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initStresses()
{
   StressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   StressYX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   StressZX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressZY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initSumOfStresses()
{
   SumStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   SumStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   SumStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   SumStressYX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   SumStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   SumStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   SumStressZX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   SumStressZY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   SumStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanOfStresses()
{
   meanStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   meanStressYX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   meanStressZX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanStressZY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initPlanarAveragingMQ()
{
   PlanarVx.resize(dimensions[2], 0);
   PlanarVy.resize(dimensions[2], 0);
   PlanarVz.resize(dimensions[2], 0);
   PlanarPr.resize(dimensions[2], 0);

   PlanarStressXX.resize(dimensions[2], 0);
   PlanarStressXY.resize(dimensions[2], 0);
   PlanarStressXZ.resize(dimensions[2], 0);

   PlanarStressYX.resize(dimensions[2], 0);
   PlanarStressYY.resize(dimensions[2], 0);
   PlanarStressYZ.resize(dimensions[2], 0);

   PlanarStressZX.resize(dimensions[2], 0);
   PlanarStressZY.resize(dimensions[2], 0);
   PlanarStressZZ.resize(dimensions[2], 0);
}
void Averaging::sumOfVolumeAveragingValues()
{
   vector<double>& vxSum = sumVaVxMatrix.getDataVector();
   vector<double>& vySum = sumVaVyMatrix.getDataVector();
   vector<double>& vzSum = sumVaVzMatrix.getDataVector();
   vector<double>& prSum = sumVaPrMatrix.getDataVector();

   vector<double>& vxVa = vaVxMatrix.getDataVector();
   vector<double>& vyVa = vaVyMatrix.getDataVector();
   vector<double>& vzVa = vaVzMatrix.getDataVector();
   vector<double>& prVa = vaPrMatrix.getDataVector();

   int size = (int)vxVa.size();

   for (int i = 0; i < size; i++)
   {

      vxSum[i] += vxVa[i];
      vySum[i] += vyVa[i];
      vzSum[i] += vzVa[i];
      prSum[i] += prVa[i];
   }
}
void Averaging::writeVolumeAveragingValuesToBinaryFiles(std::string ffname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(vaVxMatrix.getDataVector(), ffname + "Vx" + UbSystem::toString(timeStep)+".bin");
   writeMatrixToBinaryFiles<double>(vaVyMatrix.getDataVector(), ffname + "Vy" + UbSystem::toString(timeStep)+".bin");
   writeMatrixToBinaryFiles<double>(vaVzMatrix.getDataVector(), ffname + "Vz" + UbSystem::toString(timeStep)+".bin");
   writeMatrixToBinaryFiles<double>(vaPrMatrix.getDataVector(), ffname + "Pr" + UbSystem::toString(timeStep)+".bin");
}
void Averaging::meanOfVolumeAveragingValues(int numberOfTimeSteps)
{
   vector<double>& vxSum = sumVaVxMatrix.getDataVector();
   vector<double>& vySum = sumVaVyMatrix.getDataVector();
   vector<double>& vzSum = sumVaVzMatrix.getDataVector();
   vector<double>& prSum = sumVaPrMatrix.getDataVector();
                     
   vector<double>& vxMean = meanVxMatrix.getDataVector();
   vector<double>& vyMean = meanVyMatrix.getDataVector();
   vector<double>& vzMean = meanVzMatrix.getDataVector();
   vector<double>& prMean = meanPrMatrix.getDataVector();

   int size = (int)vxSum.size();

   for (int i = 0; i < size; i++)
   {
      vxMean[i] = vxSum[i] / numberOfTimeSteps;
      vyMean[i] = vySum[i] / numberOfTimeSteps;
      vzMean[i] = vzSum[i] / numberOfTimeSteps;
      prMean[i] = prSum[i] / numberOfTimeSteps;
   }
}
void Averaging::writeMeanVolumeAveragingValuesToBinaryFiles(std::string ffname)
{
   writeMatrixToBinaryFiles<double>(vaVxMatrix.getDataVector(), ffname + "Vx" + ".bin");
   writeMatrixToBinaryFiles<double>(vaVyMatrix.getDataVector(), ffname + "Vy" + ".bin");
   writeMatrixToBinaryFiles<double>(vaVzMatrix.getDataVector(), ffname + "Vz" + ".bin");
   writeMatrixToBinaryFiles<double>(vaPrMatrix.getDataVector(), ffname + "Pr" + ".bin");
}
void Averaging::fluctuationsOfVolumeAveragingValue()
{
   vector<double>& vxVa = vaVxMatrix.getDataVector();
   vector<double>& vyVa = vaVyMatrix.getDataVector();
   vector<double>& vzVa = vaVzMatrix.getDataVector();
   vector<double>& prVa = vaPrMatrix.getDataVector();

   vector<double>& vxMean = meanVxMatrix.getDataVector();
   vector<double>& vyMean = meanVyMatrix.getDataVector();
   vector<double>& vzMean = meanVzMatrix.getDataVector();
   vector<double>& prMean = meanPrMatrix.getDataVector();

   vector<double>& vxFluc = FlucVxMatrix.getDataVector();
   vector<double>& vyFluc = FlucVyMatrix.getDataVector();
   vector<double>& vzFluc = FlucVzMatrix.getDataVector();
   vector<double>& prFluc = FlucPrMatrix.getDataVector();

   vector<double>& XXStress = StressXX.getDataVector();
   vector<double>& XYStress = StressXY.getDataVector();
   vector<double>& XZStress = StressXZ.getDataVector();

   vector<double>& YXStress = StressYX.getDataVector();
   vector<double>& YYStress = StressYY.getDataVector();
   vector<double>& YZStress = StressYZ.getDataVector();

   vector<double>& ZXStress = StressZX.getDataVector();
   vector<double>& ZYStress = StressZY.getDataVector();
   vector<double>& ZZStress = StressZZ.getDataVector();


   int size = (int)vxVa.size();
   for (int i = 0; i < size; i++)
      {
         vxFluc[i] = vxVa[i] - vxMean[i];
         vyFluc[i] = vyVa[i] - vyMean[i];
         vzFluc[i] = vzVa[i] - vzMean[i];
         prFluc[i] = prVa[i] - prMean[i];

         XXStress[i] = vxFluc[i] * vxFluc[i];
         XYStress[i] = vxFluc[i] * vyFluc[i];
         XZStress[i] = vxFluc[i] * vzFluc[i];

         YXStress[i] = vyFluc[i] * vxFluc[i];
         YYStress[i] = vyFluc[i] * vyFluc[i];
         YZStress[i] = vyFluc[i] * vzFluc[i];
        
         ZXStress[i] = vzFluc[i] * vxFluc[i];
         ZYStress[i] = vzFluc[i] * vyFluc[i];
         ZZStress[i] = vzFluc[i] * vzFluc[i];
      }                       
}
void Averaging::writeFluctuationsToBinaryFiles(std::string fname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(FlucVxMatrix.getDataVector(), fname + "Vx" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(FlucVyMatrix.getDataVector(), fname + "Vy" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(FlucVzMatrix.getDataVector(), fname + "Vz" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(FlucPrMatrix.getDataVector(), fname + "Pr" + UbSystem::toString(timeStep) + ".bin");
}
void Averaging::writeStressesToBinaryFiles(std::string fname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(StressXX.getDataVector(), fname + UbSystem::toString(timeStep) + "XX" + ".bin");
   writeMatrixToBinaryFiles<double>(StressXY.getDataVector(), fname + UbSystem::toString(timeStep) + "XY" + ".bin");
   writeMatrixToBinaryFiles<double>(StressXZ.getDataVector(), fname + UbSystem::toString(timeStep) + "XZ" + ".bin");
   writeMatrixToBinaryFiles<double>(StressYX.getDataVector(), fname + UbSystem::toString(timeStep) + "YX" + ".bin");
   writeMatrixToBinaryFiles<double>(StressYY.getDataVector(), fname + UbSystem::toString(timeStep) + "YY" + ".bin");
   writeMatrixToBinaryFiles<double>(StressYZ.getDataVector(), fname + UbSystem::toString(timeStep) + "YZ" + ".bin");
   writeMatrixToBinaryFiles<double>(StressZX.getDataVector(), fname + UbSystem::toString(timeStep) + "ZX" + ".bin");
   writeMatrixToBinaryFiles<double>(StressZY.getDataVector(), fname + UbSystem::toString(timeStep) + "ZY" + ".bin");
   writeMatrixToBinaryFiles<double>(StressZZ.getDataVector(), fname + UbSystem::toString(timeStep) + "ZZ" + ".bin");
}
void Averaging::sumOfFluctuations()
{
   vector<double>& vxFluc = FlucVxMatrix.getDataVector();
   vector<double>& vyFluc = FlucVyMatrix.getDataVector();
   vector<double>& vzFluc = FlucVzMatrix.getDataVector();
   vector<double>& prFluc = FlucPrMatrix.getDataVector();

   vector<double>& SumFlucVx = sumVaVxMatrix.getDataVector();
   vector<double>& SumFlucVy = sumVaVyMatrix.getDataVector();
   vector<double>& SumFlucVz = sumVaVzMatrix.getDataVector();
   vector<double>& SumFlucPr = sumVaPrMatrix.getDataVector();

   int size = (int)vxFluc.size();

   for (int i = 0; i < size; i++)
   {

      SumFlucVx[i] += vxFluc[i];
      SumFlucVy[i] += vyFluc[i];
      SumFlucVz[i] += vzFluc[i];
      SumFlucPr[i] += prFluc[i];
   }
}

void Averaging::meanOfFluctuations(int numberOfTimeSteps)
{
   vector<double>& MeanFlucVx = meanFlucVx.getDataVector();
   vector<double>& MeanFlucVy = meanFlucVy.getDataVector();
   vector<double>& MeanFlucVz = meanFlucVz.getDataVector();
   vector<double>& MeanFlucPr = meanFlucPr.getDataVector();

   vector<double>& SumFlucVx = sumVaVxMatrix.getDataVector();
   vector<double>& SumFlucVy = sumVaVyMatrix.getDataVector();
   vector<double>& SumFlucVz = sumVaVzMatrix.getDataVector();
   vector<double>& SumFlucPr = sumVaPrMatrix.getDataVector();

   int size = (int)SumFlucVx.size();

   for (int i = 0; i < size; i++)
   {
      MeanFlucVx[i] = SumFlucVx[i] / numberOfTimeSteps;
      MeanFlucVy[i] = SumFlucVy[i] / numberOfTimeSteps;
      MeanFlucVz[i] = SumFlucVz[i] / numberOfTimeSteps;
      MeanFlucPr[i] = SumFlucPr[i] / numberOfTimeSteps;
   }
}
void Averaging::SumOfStresses()
{
   vector<double>& XXStress = StressXX.getDataVector();
   vector<double>& XYStress = StressXY.getDataVector();
   vector<double>& XZStress = StressXZ.getDataVector();

   vector<double>& YXStress = StressYX.getDataVector();
   vector<double>& YYStress = StressYY.getDataVector();
   vector<double>& YZStress = StressYZ.getDataVector();

   vector<double>& ZXStress = StressZX.getDataVector();
   vector<double>& ZYStress = StressZY.getDataVector();
   vector<double>& ZZStress = StressZZ.getDataVector();

   vector<double>& XXSum = SumStressXX.getDataVector();
   vector<double>& XYSum = SumStressXY.getDataVector();
   vector<double>& XZSum = SumStressXZ.getDataVector();
                           
   vector<double>& YXSum = SumStressYX.getDataVector();
   vector<double>& YYSum = SumStressYY.getDataVector();
   vector<double>& YZSum = SumStressYZ.getDataVector();
                           
   vector<double>& ZXSum = SumStressZX.getDataVector();
   vector<double>& ZYSum = SumStressZY.getDataVector();
   vector<double>& ZZSum = SumStressZZ.getDataVector();

   int size = (int)XXStress.size();

   for (int i = 0; i < size; i++)
   {

      XXSum[i] += XXStress[i];
      XYSum[i] += XYStress[i];
      XZSum[i] += XZStress[i];

      YXSum[i] += YXStress[i];
      YYSum[i] += YYStress[i];
      YZSum[i] += YZStress[i];

      ZXSum[i] += ZXStress[i];
      ZYSum[i] += ZYStress[i];
      ZZSum[i] += ZZStress[i];

   }
}
void Averaging::MeanOfStresses(int numberOfTimeSteps)
{
   vector<double>& XXSum = SumStressXX.getDataVector();
   vector<double>& XYSum = SumStressXY.getDataVector();
   vector<double>& XZSum = SumStressXZ.getDataVector();

   vector<double>& YXSum = SumStressYX.getDataVector();
   vector<double>& YYSum = SumStressYY.getDataVector();
   vector<double>& YZSum = SumStressYZ.getDataVector();

   vector<double>& ZXSum = SumStressZX.getDataVector();
   vector<double>& ZYSum = SumStressZY.getDataVector();
   vector<double>& ZZSum = SumStressZZ.getDataVector();

   vector<double>& XXMean = meanStressXX.getDataVector();
   vector<double>& XYMean = meanStressXY.getDataVector();
   vector<double>& XZMean = meanStressXZ.getDataVector();

   vector<double>& YXMean = meanStressYX.getDataVector();
   vector<double>& YYMean = meanStressYY.getDataVector();
   vector<double>& YZMean = meanStressYZ.getDataVector();

   vector<double>& ZXMean = meanStressZX.getDataVector();
   vector<double>& ZYMean = meanStressZY.getDataVector();
   vector<double>& ZZMean = meanStressZZ.getDataVector();

   int size = (int)XXSum.size();

   for (int i = 0; i < size; i++)
   {
      XXMean[i] = XXSum[i] / numberOfTimeSteps;
      XYMean[i] = XYSum[i] / numberOfTimeSteps;
      XZMean[i] = XZSum[i] / numberOfTimeSteps;

      XXMean[i] = XXSum[i] / numberOfTimeSteps;
      XYMean[i] = XYSum[i] / numberOfTimeSteps;
      XZMean[i] = XZSum[i] / numberOfTimeSteps;

      XXMean[i] = XXSum[i] / numberOfTimeSteps;
      XYMean[i] = XYSum[i] / numberOfTimeSteps;
      XZMean[i] = XZSum[i] / numberOfTimeSteps;
     
   }
}
void Averaging::PlanarAveragingMQ(std::array<int, 3> dimensions)
{
   double numberof_XY_points = dimensions[0] * dimensions[1];

   for (int z = 0; z < dimensions[2]; z++)
   {
      double sumVx = 0, sumVy = 0, sumVz = 0, sumPr = 0;
      double sumStressXX = 0, sumStressXY = 0, sumStressXZ = 0;
      double sumStressYX = 0, sumStressYY = 0, sumStressYZ = 0;
      double sumStressZX = 0, sumStressZY = 0, sumStressZZ = 0;
      for (int y = 0; y < dimensions[1]; y++)
         for (int x = 0; x < dimensions[0]; x++)
         {
            sumVx += meanVxMatrix(x, y, z);
            sumVy += meanVyMatrix(x, y, z);
            sumVz += meanVzMatrix(x, y, z);
            sumPr += meanPrMatrix(x, y, z);

            sumStressXX += StressXX(x, y, z);
            sumStressXY += StressXY(x, y, z);
            sumStressXZ += StressXZ(x, y, z);

            sumStressYX += StressYX(x, y, z);
            sumStressYY += StressYY(x, y, z);
            sumStressYZ += StressYZ(x, y, z);

            sumStressZX += StressZX(x, y, z);
            sumStressZY += StressZY(x, y, z);
            sumStressZZ += StressZZ(x, y, z);
         }
      PlanarVx[z] = sumVx / numberof_XY_points;
      PlanarVy[z] = sumVy / numberof_XY_points;
      PlanarVz[z] = sumVz / numberof_XY_points;
      PlanarPr[z] = sumPr / numberof_XY_points;

      PlanarStressXX[z] = sumStressXX / numberof_XY_points;
      PlanarStressXY[z] = sumStressXY / numberof_XY_points;
      PlanarStressXZ[z] = sumStressXZ / numberof_XY_points;
                                        
      PlanarStressYX[z] = sumStressYX / numberof_XY_points;
      PlanarStressYY[z] = sumStressYY / numberof_XY_points;
      PlanarStressYZ[z] = sumStressYZ / numberof_XY_points;
                                       
      PlanarStressZX[z] = sumStressZX / numberof_XY_points;
      PlanarStressZY[z] = sumStressZY / numberof_XY_points;
      PlanarStressZZ[z] = sumStressZZ / numberof_XY_points;

   }
}
   void Averaging::WriteToCSV(std::string path, double origin, double deltax)
   {
      std::ofstream ostr;
      std::string fname = path + "/av/" + "av" + ".csv";

      ostr.open(fname.c_str(), std::ios_base::out);
      if (!ostr)
      {
         ostr.clear();
         std::string path = UbSystem::getPathFromString(fname);
         if (path.size() > 0)
         {
            UbSystem::makeDirectory(path);
            ostr.open(fname.c_str(), std::ios_base::out);
         }
         if (!ostr) throw UbException(UB_EXARGS, "couldn't open file " + fname);
      }
      ostr << "z;PlanarVx;PlanarVy;PlanarVz;PlanarPr;PlanarStressXX;PlanarStressXY;PlanarStressXZ;PlanarStressYX;PlanarStressYY;PlanarStressYZ;PlanarStressZX;PlanarStressZY;PlanarStressZZ\n";
      for (int i = 0; i < dimensions[2]; i++)
      {
         double z = origin + (deltax*i);
         ostr << z << ";" << PlanarVx[i] << ";" << PlanarVy[i] << ";" << PlanarVz[i] << ";" << PlanarPr[i] << ";" << PlanarStressXX[i] << ";" << PlanarStressXY[i] << ";" << PlanarStressXZ[i] << ";" << PlanarStressYX[i] << ";" << PlanarStressYY[i] << ";" << PlanarStressYZ[i] << ";" << PlanarStressZX[i] << ";" << PlanarStressZY[i] << ";" << PlanarStressZZ[i] << "\n";
      }
      ostr.close();
   }
void Averaging::readVolumeAveragingValuesFromBinaryFiles(std::string fname, int timeStep)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + UbSystem::toString(timeStep) + ".bin", vaVxMatrix.getDataVector());
   readMatrixFromBinaryFiles<double>(fname + "Vy" + UbSystem::toString(timeStep) + ".bin", vaVyMatrix.getDataVector());
   readMatrixFromBinaryFiles<double>(fname + "Vz" + UbSystem::toString(timeStep) + ".bin", vaVzMatrix.getDataVector());
   readMatrixFromBinaryFiles<double>(fname + "Pr" + UbSystem::toString(timeStep) + ".bin", vaPrMatrix.getDataVector());
}
//////////////////////////////////////////////////////////////////////////
double Averaging::G(double x, double l)
{
   if (fabs(x) <= l)
      return l - fabs(x);
   else
      return 0.0;
}
//////////////////////////////////////////////////////////////////////////