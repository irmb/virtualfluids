#include "Averaging.h"
#include "UbLogger.h"
#include "MemoryUtil.h"
#include "UbSystem.h"

#include "ReadDataSet.h"

//#include "Postprocessing.h"
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>

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
#include <vtkXMLPUnstructuredGridReader.h>

#include <omp.h>
#include <mpi.h>

using namespace std;
void Averaging::createGeoMatrix(std::string dataNameG)
{
   UBLOG(logINFO, "createGeoMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> level_grid_timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> level_interp_timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO, "read data set from " << dataNameG << ": start");
   timer->StartTimer();

   vtkXMLPUnstructuredGridReader* reader = vtkXMLPUnstructuredGridReader::New();
   reader->SetFileName(dataNameG.c_str());
   reader->Update();

   UBLOG(logINFO, "read data set from " + dataNameG + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " << UbSystem::toString(timer->GetElapsedTime()) + " s");

   geoMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 1);

   UBLOG(logINFO, "Perform the solid nodes: start");
   level_interp_timer->StartTimer();

   vtkThreshold* thrFilter = vtkThreshold::New();
   thrFilter->SetInputData(reader->GetOutput());
   thrFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Geometry");
   thrFilter->ThresholdBetween(1, 1);
   thrFilter->Update();
   vtkUnstructuredGrid* ugrid = thrFilter->GetOutput();

   vtkPoints* points = vtkPoints::New();
   vtkCellArray* cells = vtkCellArray::New();

   int numberOfCells = ugrid->GetNumberOfCells();

   double x[3];
   array<double, 3> xMin;
   array<double, 3> xMax;
   array<int, 3> ixMin;
   array<int, 3> ixMax;
   vtkIdType idc = 0;

   for (int i = 0; i < numberOfCells; i++)
   {
      vtkIdList* plist = vtkIdList::New();
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

      getNodeIndexes(xMin, ixMin);
      getNodeIndexes(xMax, ixMax);

      for (int k = ixMin[2]; k <= ixMax[2]; k++)
      {
         for (int j = ixMin[1]; j <= ixMax[1]; j++)
         {
            for (int i = ixMin[0]; i <= ixMax[0]; i++)
            {
               if (i >= 0 && i < dimensions[0] && j >= 0 && j < dimensions[1] && k >= 0 && k < dimensions[2])
               {
                  geoMatrix(i, j, k) = 0;
               }
            }
         }
      }
      plist->Delete();
   }

   reader->Delete();
   thrFilter->Delete();
   points->Delete();
   cells->Delete();

   UBLOG(logINFO, "Perform the solid nodes: end");
   level_interp_timer->StopTimer();
   UBLOG(logINFO, "interpolation time: " << UbSystem::toString(level_interp_timer->GetElapsedTime()) << " s");

   UBLOG(logINFO, "createGeoMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO, "total time: " << UbSystem::toString(timer_total->GetElapsedTime()) << " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::writeGeoMatrixToImageFile(std::string output)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output + ".vti";

   UBLOG(logINFO, "write data set to " << vtkfilename << ": start");
   timer_write->StartTimer();

   vtkImageData* image = vtkImageData::New();

   image->SetExtent(&geo_extent[0]);
   image->SetOrigin(&geo_origin[0]);
   image->SetSpacing(&geo_spacing[0]);

   vtkIntArray* geoArray = vtkIntArray::New();
   geoArray->SetNumberOfComponents(1);
   geoArray->SetName("geo");

   int size = dimensions[0] * dimensions[1] * dimensions[2];

   geoArray->SetArray(geoMatrix.getStartAdressOfSortedArray(0, 0, 0), size, 1);
   image->GetPointData()->AddArray(geoArray);

   vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
   writer->SetInputData(image);
   writer->SetFileName(vtkfilename.c_str());
   //writer->SetDataModeToAscii();
   writer->SetDataModeToAppended();
   writer->SetCompressorTypeToZLib();
   writer->Update();

   image->Delete();
   geoArray->Delete();
   writer->Delete();

   UBLOG(logINFO, "write data set: end");
   timer_write->StopTimer();
   UBLOG(logINFO, "write data set time: " << UbSystem::toString(timer_write->GetElapsedTime()) << " s");
}
void Averaging::initMeanMqValues()
{
   meanVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::sumMqValues()
{
   vector<double>& vxSum = meanVxMatrix.getDataVector();
   vector<double>& vySum = meanVyMatrix.getDataVector();
   vector<double>& vzSum = meanVzMatrix.getDataVector();
   vector<double>& prSum = meanPrMatrix.getDataVector();

   vector<double>& vxVa = vxMatrix.getDataVector();
   vector<double>& vyVa = vyMatrix.getDataVector();
   vector<double>& vzVa = vzMatrix.getDataVector();
   vector<double>& prVa = prMatrix.getDataVector();

   int size = (int)vxVa.size();

   for (int i = 0; i < size; i++)
   {
      vxSum[i] += vxVa[i];
      vySum[i] += vyVa[i];
      vzSum[i] += vzVa[i];
      prSum[i] += prVa[i];
   }
}
void Averaging::computeMeanMqValues(int numberOfTimeSteps)
{
   vector<double>& vxSum = meanVxMatrix.getDataVector();
   vector<double>& vySum = meanVyMatrix.getDataVector();
   vector<double>& vzSum = meanVzMatrix.getDataVector();
   vector<double>& prSum = meanPrMatrix.getDataVector();

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
void Averaging::writeMeanMqValuesToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<double>(meanVxMatrix, fname + "Vx" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVyMatrix, fname + "Vy" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVzMatrix, fname + "Vz" + ".bin");
   writeMatrixToBinaryFiles<double>(meanPrMatrix, fname + "Pr" + ".bin");
}
void Averaging::readMeanMqValuesFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + ".bin", meanVxMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + ".bin", meanVyMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + ".bin", meanVzMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + ".bin", meanPrMatrix);
}
void Averaging::volumeAveragingOfMeanMqValuesWithMPI(double l_real)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "volume averaging: start");
   timer_averaging->StartTimer();

   double l = round(l_real / deltax);
   UBLOG(logINFO, "l = " + UbSystem::toString(l));

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) << " x " + UbSystem::toString(dimensions[1]) << " x " << UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   vaMeanVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

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

#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               int ID = omp_get_thread_num();
               if (i == 0 && ID == 0)
               {
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
               }
               double vx = 0.0;
               double vy = 0.0;
               double vz = 0.0;
               double pr = 0.0;

               int ll = (int)l;
               int llz1 = ll;
               if (x3 - ll < 0) llz1 = x3;
               if (x3 + ll >= dimensions[2]) llz1 = dimensions[2] - 1 - x3;
               double lQuadrat = l * l;
               double lNorm = lQuadrat * lQuadrat * double(llz1) * double(llz1);

               //#pragma omp parallel for
               for (int z = -llz1; z <= +llz1; z++)
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

                        if (zz < 0) zz = 0;
                        if (zz >= dimensions[2]) zz = dimensions[2] - 1;

                        double mm = (G((double)x, l) * G((double)y, l) * G((double)z, (double)llz1)) / lNorm;
                        double gamma = (double)geoMatrix(xx, yy, zz);

                        vx += gamma * mm * meanVxMatrix(xx, yy, zz);
                        vy += gamma * mm * meanVyMatrix(xx, yy, zz);
                        vz += gamma * mm * meanVzMatrix(xx, yy, zz);
                        pr += gamma * mm * meanPrMatrix(xx, yy, zz);

                     }

               vaMeanVxMatrix(x1, x2, x3) = vx;
               vaMeanVyMatrix(x1, x2, x3) = vy;
               vaMeanVzMatrix(x1, x2, x3) = vz;
               vaMeanPrMatrix(x1, x2, x3) = pr;

               if (i % p == 0 && i != 0 && ID == 0)
               {
                  timer_inloop->StopTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
                  UBLOG(logINFO, "time per " + UbSystem::toString(p) + " points: " + UbSystem::toString(timer_inloop->GetElapsedTime()) + " s");
                  UBLOG(logINFO, "actual memory usage: " << UbSystem::toString(Utilities::getPhysMemUsedByMe() / 1e9) << " GByte");
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "thread id: " + UbSystem::toString(ID));
                  UBLOG(logINFO, "Number of treads: " + UbSystem::toString(omp_get_num_threads()));
               }
               i++;
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
         for (int x3 = 0; x3 < dimensions[2]; x3++)
            for (int x2 = 0; x2 < dimensions[1]; x2++)
               for (int x1 = lstartX1; x1 < lstopX1; x1++)
               {
                  vaMeanVxMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaMeanVyMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaMeanVzMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaMeanPrMatrix(x1, x2, x3) = receiveBuffer[c++];
               }
      }
   }
   else
   {
      vector<double> sendBuffer;
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               sendBuffer.push_back(vaMeanVxMatrix(x1, x2, x3));
               sendBuffer.push_back(vaMeanVyMatrix(x1, x2, x3));
               sendBuffer.push_back(vaMeanVzMatrix(x1, x2, x3));
               sendBuffer.push_back(vaMeanPrMatrix(x1, x2, x3));
            }
      int count = (int)sendBuffer.size();
      MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   }

   timer_averaging->StopTimer();
   UBLOG(logINFO, "volume averaging: end");
   UBLOG(logINFO, "volume averaging time: " + UbSystem::toString(timer_averaging->GetElapsedTime()) + " s");
}
void Averaging::writeVaMeanMqValuesToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<double>(vaMeanVxMatrix, fname + "Vx" + ".bin");
   writeMatrixToBinaryFiles<double>(vaMeanVyMatrix, fname + "Vy" + ".bin");
   writeMatrixToBinaryFiles<double>(vaMeanVzMatrix, fname + "Vz" + ".bin");
   writeMatrixToBinaryFiles<double>(vaMeanPrMatrix, fname + "Pr" + ".bin");
}
void Averaging::readVaMeanMqValuesFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + ".bin", vaMeanVxMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + ".bin", vaMeanVyMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + ".bin", vaMeanVzMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + ".bin", vaMeanPrMatrix);
}
void Averaging::initFluctuationsOfMqValues()
{
   flucVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   flucVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   flucVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   flucPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::computeFluctuationsOfMqValues()
{
   vector<double>& vxF = vxMatrix.getDataVector();
   vector<double>& vyF = vyMatrix.getDataVector();
   vector<double>& vzF = vzMatrix.getDataVector();
   vector<double>& prF = prMatrix.getDataVector();

   vector<double>& vxMean = meanVxMatrix.getDataVector();
   vector<double>& vyMean = meanVyMatrix.getDataVector();
   vector<double>& vzMean = meanVzMatrix.getDataVector();
   vector<double>& prMean = meanPrMatrix.getDataVector();

   vector<double>& vxFluc = flucVxMatrix.getDataVector();
   vector<double>& vyFluc = flucVyMatrix.getDataVector();
   vector<double>& vzFluc = flucVzMatrix.getDataVector();
   vector<double>& prFluc = flucPrMatrix.getDataVector();

   int size = (int)vxF.size();

   for (int i = 0; i < size; i++)
   {
      vxFluc[i] = vxF[i] - vxMean[i];
      vyFluc[i] = vyF[i] - vyMean[i];
      vzFluc[i] = vzF[i] - vzMean[i];
      prFluc[i] = prF[i] - prMean[i];
   }
}
void Averaging::writeFluctuationsOfMqValuesToBinaryFiles(std::string fname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(flucVxMatrix, fname + "Vx" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(flucVyMatrix, fname + "Vy" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(flucVzMatrix, fname + "Vz" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(flucPrMatrix, fname + "Pr" + UbSystem::toString(timeStep) + ".bin");
}
void Averaging::readFluctuationsOfMqValuesFromBinaryFiles(std::string fname, int timeStep)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + UbSystem::toString(timeStep) + ".bin", flucVxMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + UbSystem::toString(timeStep) + ".bin", flucVyMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + UbSystem::toString(timeStep) + ".bin", flucVzMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + UbSystem::toString(timeStep) + ".bin", flucPrMatrix);
}
void Averaging::volumeAveragingOfFluctuationsWithMPI(double l_real)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "volume averaging fluct and stress: start");
   timer_averaging->StartTimer();

   double l = round(l_real / deltax);
   UBLOG(logINFO, "l = " + UbSystem::toString(l));

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) << " x " + UbSystem::toString(dimensions[1]) << " x " << UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   vaFlucVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

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
#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               int ID = omp_get_thread_num();
               if (i == 0 && ID == 0)
               {
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
               }
               double flucvx = 0.0;
               double flucvy = 0.0;
               double flucvz = 0.0;
               double flucpr = 0.0;

               int ll = (int)l;
               int llz1 = ll;
               if (x3 - ll < 0) llz1 = x3;
               if (x3 + ll >= dimensions[2]) llz1 = dimensions[2] - 1 - x3;
               double lQuadrat = l * l;
               double lNorm = lQuadrat * lQuadrat * double(llz1) * double(llz1);

               //#pragma omp parallel for
               for (int z = -llz1; z <= +llz1; z++)
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

                        double mm = (G((double)x, l) * G((double)y, l) * G((double)z, (double)llz1)) / lNorm;
                        double gamma = (double)geoMatrix(xx, yy, zz);

                        flucvx += gamma * mm * flucVxMatrix(xx, yy, zz);
                        flucvy += gamma * mm * flucVyMatrix(xx, yy, zz);
                        flucvz += gamma * mm * flucVzMatrix(xx, yy, zz);
                        flucpr += gamma * mm * flucPrMatrix(xx, yy, zz);
                     }

               vaFlucVxMatrix(x1, x2, x3) = flucvx;
               vaFlucVyMatrix(x1, x2, x3) = flucvy;
               vaFlucVzMatrix(x1, x2, x3) = flucvz;
               vaFlucPrMatrix(x1, x2, x3) = flucpr;

               if (i % p == 0 && i != 0 && ID == 0)
               {
                  timer_inloop->StopTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
                  UBLOG(logINFO, "time per " + UbSystem::toString(p) + " points: " + UbSystem::toString(timer_inloop->GetElapsedTime()) + " s");
                  UBLOG(logINFO, "actual memory usage: " << UbSystem::toString(Utilities::getPhysMemUsedByMe() / 1e9) << " GByte");
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "thread id: " + UbSystem::toString(ID));
                  UBLOG(logINFO, "Number of treads: " + UbSystem::toString(omp_get_num_threads()));
               }
               i++;
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
         for (int x3 = 0; x3 < dimensions[2]; x3++)
            for (int x2 = 0; x2 < dimensions[1]; x2++)
               for (int x1 = lstartX1; x1 < lstopX1; x1++)
               {
                  vaFlucVxMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaFlucVyMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaFlucVzMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaFlucPrMatrix(x1, x2, x3) = receiveBuffer[c++];
               }
      }
   }
   else
   {
      vector<double> sendBuffer;
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               sendBuffer.push_back(vaFlucVxMatrix(x1, x2, x3));
               sendBuffer.push_back(vaFlucVyMatrix(x1, x2, x3));
               sendBuffer.push_back(vaFlucVzMatrix(x1, x2, x3));
               sendBuffer.push_back(vaFlucPrMatrix(x1, x2, x3));
            }
      int count = (int)sendBuffer.size();
      MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   }

   timer_averaging->StopTimer();
   UBLOG(logINFO, "volume averaging fluct and stress: end");
   UBLOG(logINFO, "volume averaging fluct and stress time: " + UbSystem::toString(timer_averaging->GetElapsedTime()) + " s");
}
void Averaging::readTimeAveragedDataFromVtkFile(std::string dataNameMQ)
{
   UBLOG(logINFO, "createMQMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO, "read data set from " + dataNameMQ + ": start");
   timer->StartTimer();

   vtkXMLPUnstructuredGridReader* reader = vtkXMLPUnstructuredGridReader::New();
   reader->SetFileName(dataNameMQ.c_str());
   reader->Update();

   UBLOG(logINFO, "read data set from " + dataNameMQ + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " + UbSystem::toString(timer->GetElapsedTime()) + " s");

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) + " x " + UbSystem::toString(dimensions[1]) + " x " + UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   meanVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   StressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   vtkUnstructuredGrid* ugrid = reader->GetOutput();

   vtkPoints* points = vtkPoints::New();
   vtkCellArray* cells = vtkCellArray::New();

   int numberOfCells = ugrid->GetNumberOfCells();

   vtkDataArray* vxArray = ugrid->GetPointData()->GetArray("taVx");
   vtkDataArray* vyArray = ugrid->GetPointData()->GetArray("taVy");
   vtkDataArray* vzArray = ugrid->GetPointData()->GetArray("taVz");
   vtkDataArray* prArray = ugrid->GetPointData()->GetArray("taRho");

   vtkDataArray* vxxArray = ugrid->GetPointData()->GetArray("taVxx");
   vtkDataArray* vyyArray = ugrid->GetPointData()->GetArray("taVyy");
   vtkDataArray* vzzArray = ugrid->GetPointData()->GetArray("taVzz");
   vtkDataArray* vxyArray = ugrid->GetPointData()->GetArray("taVxy");
   vtkDataArray* vxzArray = ugrid->GetPointData()->GetArray("taVxz");
   vtkDataArray* vyzArray = ugrid->GetPointData()->GetArray("taVyz");
   //Vxx; Vyy; Vzz; Vxy; Vxz; Vyz

   double x[3];
   array<double, 3> xMin;
   array<double, 3> xMax;
   array<int, 3> ixMin;
   array<int, 3> ixMax;
   vtkIdType idc = 0;

   for (int i = 0; i < numberOfCells; i++)
   {
      vtkIdList* plist = vtkIdList::New();
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

      getNodeIndexes(xMin, ixMin);
      getNodeIndexes(xMax, ixMax);
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
                     meanVxMatrix(i, j, k) = vxArray->GetTuple1(plist->GetId(c));
                     meanVyMatrix(i, j, k) = vyArray->GetTuple1(plist->GetId(c));
                     meanVzMatrix(i, j, k) = vzArray->GetTuple1(plist->GetId(c));
                     meanPrMatrix(i, j, k) = prArray->GetTuple1(plist->GetId(c));

                     StressXX(i, j, k) = vxxArray->GetTuple1(plist->GetId(c));
                     StressYY(i, j, k) = vyyArray->GetTuple1(plist->GetId(c));
                     StressZZ(i, j, k) = vzzArray->GetTuple1(plist->GetId(c));
                     StressXY(i, j, k) = vxyArray->GetTuple1(plist->GetId(c));
                     StressXZ(i, j, k) = vxzArray->GetTuple1(plist->GetId(c));
                     StressYZ(i, j, k) = vyzArray->GetTuple1(plist->GetId(c));
                     c++;
                  }
               }
            }
         }
      }
      plist->Delete();
   }

   reader->Delete();
   points->Delete();
   cells->Delete();

   UBLOG(logINFO, "createMQMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO, "total time: " + UbSystem::toString(timer_total->GetElapsedTime()) + " s");
}
void Averaging::volumeAveragingOfTimeAveragedDataWithMPI(double l_real)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "volume averaging: start");
   timer_averaging->StartTimer();

   double l = round(l_real / deltax);
   UBLOG(logINFO, "l = " + UbSystem::toString(l));

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) << " x " + UbSystem::toString(dimensions[1]) << " x " << UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   vaMeanVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   vaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

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

#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               int ID = omp_get_thread_num();
               if (i == 0 && ID == 0)
               {
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
               }
               double vx = 0.0;
               double vy = 0.0;
               double vz = 0.0;
               double pr = 0.0;
               double stressXX = 0.0;
               double stressYY = 0.0;
               double stressZZ = 0.0;
               double stressXY = 0.0;
               double stressXZ = 0.0;
               double stressYZ = 0.0;

               int ll = (int)l;
               int llz1 = ll;
               if (x3 - ll < 0) llz1 = x3;
               if (x3 + ll >= dimensions[2]) llz1 = dimensions[2] - 1 - x3;
               double lQuadrat = l * l;
               double lNorm = lQuadrat * lQuadrat * double(llz1)*double(llz1);

               //#pragma omp parallel for
               for (int z = -llz1; z <= +llz1; z++)
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

                        if (zz < 0) zz = 0;
                        if (zz >= dimensions[2]) zz = dimensions[2] - 1;

                        if (x3 != 0 || x3 != dimensions[2] - 1)
                        {
                           double mm = (G((double)x, l) * G((double)y, l) * G((double)z, (double)llz1)) / lNorm;
                           double gamma = (double)geoMatrix(xx, yy, zz);

                           vx += gamma * mm * meanVxMatrix(xx, yy, zz);
                           vy += gamma * mm * meanVyMatrix(xx, yy, zz);
                           vz += gamma * mm * meanVzMatrix(xx, yy, zz);
                           pr += gamma * mm * meanPrMatrix(xx, yy, zz);

                           stressXX += gamma * mm * StressXX(xx, yy, zz);
                           stressYY += gamma * mm * StressYY(xx, yy, zz);
                           stressZZ += gamma * mm * StressZZ(xx, yy, zz);
                           stressXY += gamma * mm * StressXY(xx, yy, zz);
                           stressXZ += gamma * mm * StressXZ(xx, yy, zz);
                           stressYZ += gamma * mm * StressYZ(xx, yy, zz);
                        }
                        else
                        {
                           vx += meanVxMatrix(xx, yy, zz);
                           vy += meanVyMatrix(xx, yy, zz);
                           vz += meanVzMatrix(xx, yy, zz);
                           pr += meanPrMatrix(xx, yy, zz);

                           stressXX += StressXX(xx, yy, zz);
                           stressYY += StressYY(xx, yy, zz);
                           stressZZ += StressZZ(xx, yy, zz);
                           stressXY += StressXY(xx, yy, zz);
                           stressXZ += StressXZ(xx, yy, zz);
                           stressYZ += StressYZ(xx, yy, zz);
                        }
                     }

               vaMeanVxMatrix(x1, x2, x3) = vx;
               vaMeanVyMatrix(x1, x2, x3) = vy;
               vaMeanVzMatrix(x1, x2, x3) = vz;
               vaMeanPrMatrix(x1, x2, x3) = pr;

               vaStressXX(x1, x2, x3) = stressXX;
               vaStressYY(x1, x2, x3) = stressYY;
               vaStressZZ(x1, x2, x3) = stressZZ;
               vaStressXY(x1, x2, x3) = stressXY;
               vaStressXZ(x1, x2, x3) = stressXZ;
               vaStressYZ(x1, x2, x3) = stressYZ;

               if (i % p == 0 && i != 0 && ID == 0)
               {
                  timer_inloop->StopTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
                  UBLOG(logINFO, "time per " + UbSystem::toString(p) + " points: " + UbSystem::toString(timer_inloop->GetElapsedTime()) + " s");
                  UBLOG(logINFO, "actual memory usage: " << UbSystem::toString(Utilities::getPhysMemUsedByMe() / 1e9) << " GByte");
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "thread id: " + UbSystem::toString(ID));
                  UBLOG(logINFO, "Number of treads: " + UbSystem::toString(omp_get_num_threads()));
               }
               i++;
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
         for (int x3 = 0; x3 < dimensions[2]; x3++)
            for (int x2 = 0; x2 < dimensions[1]; x2++)
               for (int x1 = lstartX1; x1 < lstopX1; x1++)
               {
                  vaMeanVxMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaMeanVyMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaMeanVzMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaMeanPrMatrix(x1, x2, x3) = receiveBuffer[c++];

                  vaStressXX(x1, x2, x3) = receiveBuffer[c++];
                  vaStressYY(x1, x2, x3) = receiveBuffer[c++];
                  vaStressZZ(x1, x2, x3) = receiveBuffer[c++];
                  vaStressXY(x1, x2, x3) = receiveBuffer[c++];
                  vaStressXZ(x1, x2, x3) = receiveBuffer[c++];
                  vaStressYZ(x1, x2, x3) = receiveBuffer[c++];
               }
      }
   }
   else
   {
      vector<double> sendBuffer;
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               sendBuffer.push_back(vaMeanVxMatrix(x1, x2, x3));
               sendBuffer.push_back(vaMeanVyMatrix(x1, x2, x3));
               sendBuffer.push_back(vaMeanVzMatrix(x1, x2, x3));
               sendBuffer.push_back(vaMeanPrMatrix(x1, x2, x3));

               sendBuffer.push_back(vaStressXX(x1, x2, x3));
               sendBuffer.push_back(vaStressYY(x1, x2, x3));
               sendBuffer.push_back(vaStressZZ(x1, x2, x3));
               sendBuffer.push_back(vaStressXY(x1, x2, x3));
               sendBuffer.push_back(vaStressXZ(x1, x2, x3));
               sendBuffer.push_back(vaStressYZ(x1, x2, x3));
            }
      int count = (int)sendBuffer.size();
      MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   }

   timer_averaging->StopTimer();
   UBLOG(logINFO, "volume averaging: end");
   UBLOG(logINFO, "volume averaging time: " + UbSystem::toString(timer_averaging->GetElapsedTime()) + " s");
}
void Averaging::planarAveragingOfVaTaData()
{
   double numberof_XY_points = (double)dimensions[0] * (double)dimensions[1];

   for (int z = 0; z < dimensions[2]; z++)
   {
      double sumVx = 0, sumVy = 0, sumVz = 0, sumPr = 0;
      double sumFluctVx = 0, sumFluctVy = 0, sumFluctVz = 0, sumFluctPr = 0;
      double sumStressXX = 0, sumStressYY = 0, sumStressZZ = 0, sumStressXY = 0, sumStressXZ = 0, sumStressYZ = 0;
      for (int y = 0; y < dimensions[1]; y++)
         for (int x = 0; x < dimensions[0]; x++)
         {
            sumVx += vaMeanVxMatrix(x, y, z);
            sumVy += vaMeanVyMatrix(x, y, z);
            sumVz += vaMeanVzMatrix(x, y, z);
            sumPr += vaMeanPrMatrix(x, y, z);

            sumStressXX += vaStressXX(x, y, z);
            sumStressYY += vaStressYY(x, y, z);
            sumStressZZ += vaStressZZ(x, y, z);
            sumStressXY += vaStressXY(x, y, z);
            sumStressXZ += vaStressXZ(x, y, z);
            sumStressYZ += vaStressYZ(x, y, z);
         }
      PlanarVx[z] = sumVx / numberof_XY_points;
      PlanarVy[z] = sumVy / numberof_XY_points;
      PlanarVz[z] = sumVz / numberof_XY_points;
      PlanarPr[z] = sumPr / numberof_XY_points;

      PlanarStressXX[z] = sumStressXX / numberof_XY_points;
      PlanarStressYY[z] = sumStressYY / numberof_XY_points;
      PlanarStressZZ[z] = sumStressZZ / numberof_XY_points;
      PlanarStressXY[z] = sumStressXY / numberof_XY_points;
      PlanarStressXZ[z] = sumStressXZ / numberof_XY_points;
      PlanarStressYZ[z] = sumStressYZ / numberof_XY_points;
   }
}
//////////////////////////////////////////////////////////////////////////
void Averaging::getNodeIndexes(std::array<double, 3> x, std::array<int, 3>& ix)
{
   ix[0] = (int)round((x[0] - geo_origin[0]) / deltax);
   ix[1] = (int)round((x[1] - geo_origin[1]) / deltax);
   ix[2] = (int)round((x[2] - geo_origin[2]) / deltax);
}
//////////////////////////////////////////////////////////////////////////
void Averaging::createMQMatrix(std::string dataNameMQ)
{
   UBLOG(logINFO, "createMQMatrix:start");

   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
   vtkSmartPointer<vtkTimerLog> timer_total = vtkSmartPointer<vtkTimerLog>::New();

   timer_total->StartTimer();

   UBLOG(logINFO, "read data set from " + dataNameMQ + ": start");
   timer->StartTimer();

   vtkXMLPUnstructuredGridReader* reader = vtkXMLPUnstructuredGridReader::New();
   reader->SetFileName(dataNameMQ.c_str());
   reader->Update();

   UBLOG(logINFO, "read data set from " + dataNameMQ + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " + UbSystem::toString(timer->GetElapsedTime()) + " s");

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) + " x " + UbSystem::toString(dimensions[1]) + " x " + UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   vxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   prMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   vtkUnstructuredGrid* ugrid = reader->GetOutput(); 

   vtkPoints* points = vtkPoints::New();
   vtkCellArray* cells = vtkCellArray::New();

   int numberOfCells = ugrid->GetNumberOfCells();

   vtkDataArray* vxArray = ugrid->GetPointData()->GetArray("Vx");
   vtkDataArray* vyArray = ugrid->GetPointData()->GetArray("Vy");
   vtkDataArray* vzArray = ugrid->GetPointData()->GetArray("Vz");
   vtkDataArray* prArray = ugrid->GetPointData()->GetArray("Rho");

   double x[3];
   array<double, 3> xMin;
   array<double, 3> xMax;
   array<int, 3> ixMin;
   array<int, 3> ixMax;
   vtkIdType idc = 0;

   for (int i = 0; i < numberOfCells; i++)
   {
      vtkIdList* plist = vtkIdList::New();
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

      getNodeIndexes(xMin, ixMin);
      getNodeIndexes(xMax, ixMax);
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
      plist->Delete();
   }

   reader->Delete();
   points->Delete();
   cells->Delete();

   UBLOG(logINFO, "createMQMatrix:end");
   timer_total->StopTimer();
   UBLOG(logINFO, "total time: " + UbSystem::toString(timer_total->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
void Averaging::writeMatrixToImageFile(std::string output, std::array<CbArray3D<double>, 4> matrix)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   std::string vtkfilename = output + ".vti";

   UBLOG(logINFO, "write data set to " + vtkfilename + ": start");
   timer_write->StartTimer();

   vtkImageData* image = vtkImageData::New();

   image->SetExtent(&geo_extent[0]);
   image->SetOrigin(&geo_origin[0]);
   image->SetSpacing(&geo_spacing[0]);

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

   vxArray->SetArray(matrix[0].getStartAdressOfSortedArray(0, 0, 0), size, 1);
   vyArray->SetArray(matrix[1].getStartAdressOfSortedArray(0, 0, 0), size, 1);
   vzArray->SetArray(matrix[2].getStartAdressOfSortedArray(0, 0, 0), size, 1);
   prArray->SetArray(matrix[3].getStartAdressOfSortedArray(0, 0, 0), size, 1);

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
   
   image->Delete();
   vxArray->Delete();
   vyArray->Delete();
   vzArray->Delete();
   writer->Delete();

   UBLOG(logINFO, "write data set: end");
   timer_write->StopTimer();
   UBLOG(logINFO, "write data set time: " + UbSystem::toString(timer_write->GetElapsedTime()) + " s");
}
void Averaging::writeMqMatrixToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { vxMatrix, vyMatrix, vzMatrix, prMatrix };
   writeMatrixToImageFile(output, matrix);
}
void Averaging::writeVaMatrixToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { vaVxMatrix, vaVyMatrix, vaVzMatrix, vaPrMatrix };
   writeMatrixToImageFile(output, matrix);
}
void Averaging::writeVaSumMatrixToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { sumVaVxMatrix, sumVaVyMatrix, sumVaVzMatrix, sumVaPrMatrix };
   writeMatrixToImageFile(output, matrix);
}
void Averaging::writeMeanMatrixToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { vaMeanVxMatrix, vaMeanVyMatrix, vaMeanVzMatrix, vaMeanPrMatrix };
   writeMatrixToImageFile(output, matrix);
}

//////////////////////////////////////////////////////////////////////////
void Averaging::readGeoMatrix(string dataNameG)
{
   vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "readGeoMatrix:start");

   UBLOG(logINFO, "read data set from " + dataNameG + ": start");
   timer->StartTimer();
   vtkDataSet* dataSetGeo(ReadDataSet(dataNameG.c_str()));
   UBLOG(logINFO, "read data set from " + dataNameG + ": end");
   timer->StopTimer();
   UBLOG(logINFO, "read data set time: " + UbSystem::toString(timer->GetElapsedTime()) + " s");

   vtkImageData* image = vtkImageData::SafeDownCast(dataSetGeo);

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

   vtkDataArray* geoArray = dataSetGeo->GetPointData()->GetArray("geo");

   int numberOfPoints = dataSetGeo->GetNumberOfPoints();
   int* gm = geoMatrix.getStartAdressOfSortedArray(0, 0, 0);
   for (int i = 0; i < numberOfPoints; i++)
   {
      gm[i] = (int)geoArray->GetTuple1(i);
   }

   dataSetGeo->Delete();
   image->Delete();
   geoArray->Delete();

   UBLOG(logINFO, "readGeoMatrix:end");
}
void Averaging::writeGeoMatrixToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<int>(geoMatrix, fname);
}
void Averaging::readGeoMatrixFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<int>(fname, geoMatrix);
}
void Averaging::writeMqMatrixToBinaryFiles(std::string fname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(vxMatrix, fname + "Vx" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vyMatrix, fname + "Vy" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vzMatrix, fname + "Vz" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(prMatrix, fname + "Pr" + UbSystem::toString(timeStep) + ".bin");
}
void Averaging::readMqMatrixFromBinaryFiles(std::string fname, int timeStep)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + UbSystem::toString(timeStep) + ".bin", vxMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + UbSystem::toString(timeStep) + ".bin", vyMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + UbSystem::toString(timeStep) + ".bin", vzMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + UbSystem::toString(timeStep) + ".bin", prMatrix);
}

//-------------------------------- volume avaraging --------------------------
void Averaging::initVolumeAveragingValues()
{
   sumVaVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initVolumeAveragingFluctStressValues()
{
   sumVaFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   
   sumVaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanVolumeAveragingValues()
{
   vaMeanVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaMeanPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanVolumeAveragingFluctStressValues()
{
   meanVaFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   
   meanVaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
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
void Averaging::meanOfVolumeAveragingValues(int numberOfTimeSteps)
{
   vector<double>& vxSum = sumVaVxMatrix.getDataVector();
   vector<double>& vySum = sumVaVyMatrix.getDataVector();
   vector<double>& vzSum = sumVaVzMatrix.getDataVector();
   vector<double>& prSum = sumVaPrMatrix.getDataVector();

   vector<double>& vxMean = vaMeanVxMatrix.getDataVector();
   vector<double>& vyMean = vaMeanVyMatrix.getDataVector();
   vector<double>& vzMean = vaMeanVzMatrix.getDataVector();
   vector<double>& prMean = vaMeanPrMatrix.getDataVector();

   int size = (int)vxSum.size();

   for (int i = 0; i < size; i++)
   {
      vxMean[i] = vxSum[i] / numberOfTimeSteps;
      vyMean[i] = vySum[i] / numberOfTimeSteps;
      vzMean[i] = vzSum[i] / numberOfTimeSteps;
      prMean[i] = prSum[i] / numberOfTimeSteps;
   }
}
void Averaging::volumeAveragingWithMPI(double l_real)
{
   //////////////////////////////////////////////////////////////////////////
   //DEBUG
   //////////////////////////////////////////////////////////////////////////
   //vaVxMatrix = vxMatrix;
   //vaVyMatrix = vyMatrix;
   //vaVzMatrix = vzMatrix;
   //vaPrMatrix = prMatrix;
   //return;
   //////////////////////////////////////////////////////////////////////////

   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "volume averaging: start");
   timer_averaging->StartTimer();

   double l = round(l_real / deltax);
   UBLOG(logINFO, "l = " + UbSystem::toString(l));

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
      
#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               int ID = omp_get_thread_num();
               if (i == 0 && ID == 0)
               {
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
               }
               double vx = 0.0;
               double vy = 0.0;
               double vz = 0.0;
               double pr = 0.0;

               int ll = (int)l;
               int llz1 = ll;
               if (x3 - ll < 0) llz1 = x3;
               if (x3 + ll >= dimensions[2]) llz1 = dimensions[2] - 1 - x3;
               double lQuadrat = l * l;
               double lNorm = lQuadrat * lQuadrat * double(llz1) * double(llz1);

               //#pragma omp parallel for
               for (int z = -llz1; z <= +llz1; z++)
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

                        if (zz < 0) zz = 0;
                        if (zz >= dimensions[2]) zz = dimensions[2] - 1;

                        double mm = (G((double)x, l)*G((double)y, l)*G((double)z, (double)llz1)) / lNorm;
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

               if (i%p == 0 && i != 0 && ID == 0)
               {
                  timer_inloop->StopTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
                  UBLOG(logINFO, "time per " + UbSystem::toString(p) + " points: " + UbSystem::toString(timer_inloop->GetElapsedTime()) + " s");
                  UBLOG(logINFO, "actual memory usage: " << UbSystem::toString(Utilities::getPhysMemUsedByMe() / 1e9) << " GByte");
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "thread id: " + UbSystem::toString(ID));
                  UBLOG(logINFO, "Number of treads: " + UbSystem::toString(omp_get_num_threads()));
               }
               i++;
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
         for (int x3 = 0; x3 < dimensions[2]; x3++)
            for (int x2 = 0; x2 < dimensions[1]; x2++)
               for (int x1 = lstartX1; x1 < lstopX1; x1++)
               {
                  vaVxMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaVyMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaVzMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaPrMatrix(x1, x2, x3) = receiveBuffer[c++];
               }
      }
   }
   else
   {
      vector<double> sendBuffer;
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               sendBuffer.push_back(vaVxMatrix(x1, x2, x3));
               sendBuffer.push_back(vaVyMatrix(x1, x2, x3));
               sendBuffer.push_back(vaVzMatrix(x1, x2, x3));
               sendBuffer.push_back(vaPrMatrix(x1, x2, x3));
            }
      int count = (int)sendBuffer.size();
      MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   }

   timer_averaging->StopTimer();
   UBLOG(logINFO, "volume averaging: end");
   UBLOG(logINFO, "volume averaging time: " + UbSystem::toString(timer_averaging->GetElapsedTime()) + " s");
}
void Averaging::volumeAveragingFluctStressWithMPI(double l_real)
{
   vtkSmartPointer<vtkTimerLog> timer_averaging = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO, "volume averaging fluct and stress: start");
   timer_averaging->StartTimer();

   double l = round(l_real / deltax);
   UBLOG(logINFO, "l = " + UbSystem::toString(l));

   UBLOG(logINFO, "NX1 x NX2 x NX3 = " + UbSystem::toString(dimensions[0]) << " x " + UbSystem::toString(dimensions[1]) << " x " << UbSystem::toString(dimensions[2]));

   int size = dimensions[0] * dimensions[1] * dimensions[2];
   vaFlucVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);

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
#pragma omp parallel for //private(i)//scheduler(dynamic, 1)
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               int ID = omp_get_thread_num();
               if (i == 0 && ID == 0)
               {
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
               }
               double flucvx = 0.0;
               double flucvy = 0.0;
               double flucvz = 0.0;
               double flucpr = 0.0;
               double stressXX = 0.0;
               double stressYY = 0.0;
               double stressZZ = 0.0;
               double stressXY = 0.0;
               double stressXZ = 0.0;
               double stressYZ = 0.0;

               int ll = (int)l;
               int llz1 = ll;
               if (x3 - ll < 0) llz1 = x3;
               if (x3 + ll >= dimensions[2]) llz1 = dimensions[2] - 1 - x3;
               double lQuadrat = l * l;
               double lNorm = lQuadrat * lQuadrat * double(llz1) * double(llz1);

               //#pragma omp parallel for
               for (int z = -llz1; z <= +llz1; z++)
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

                        double mm = (G((double)x, l)*G((double)y, l)*G((double)z, (double)llz1)) / lNorm;
                        double gamma = (double)geoMatrix(xx, yy, zz);

                        flucvx += gamma*mm*flucVxMatrix(xx, yy, zz);
                        flucvy += gamma*mm*flucVyMatrix(xx, yy, zz);
                        flucvz += gamma*mm*flucVzMatrix(xx, yy, zz);
                        flucpr += gamma*mm*flucPrMatrix(xx, yy, zz);

                        stressXX += gamma*mm*StressXX(xx, yy, zz);
                        stressYY += gamma*mm*StressYY(xx, yy, zz);
                        stressZZ += gamma*mm*StressZZ(xx, yy, zz);
                        stressXY += gamma*mm*StressXY(xx, yy, zz);
                        stressXZ += gamma*mm*StressXZ(xx, yy, zz);
                        stressYZ += gamma*mm*StressYZ(xx, yy, zz);

                     }

               vaFlucVxMatrix(x1, x2, x3) = flucvx;
               vaFlucVyMatrix(x1, x2, x3) = flucvy;
               vaFlucVzMatrix(x1, x2, x3) = flucvz;
               vaFlucPrMatrix(x1, x2, x3) = flucpr;

               vaStressXX(x1, x2, x3) = stressXX;
               vaStressYY(x1, x2, x3) = stressYY;
               vaStressZZ(x1, x2, x3) = stressZZ;
               vaStressXY(x1, x2, x3) = stressXY;
               vaStressXZ(x1, x2, x3) = stressXZ;
               vaStressYZ(x1, x2, x3) = stressYZ;

               if (i%p == 0 && i != 0 && ID == 0)
               {
                  timer_inloop->StopTimer();
                  UBLOG(logINFO, "point id = " + UbSystem::toString(i));
                  UBLOG(logINFO, "time per " + UbSystem::toString(p) + " points: " + UbSystem::toString(timer_inloop->GetElapsedTime()) + " s");
                  UBLOG(logINFO, "actual memory usage: " << UbSystem::toString(Utilities::getPhysMemUsedByMe() / 1e9) << " GByte");
                  timer_inloop->StartTimer();
                  UBLOG(logINFO, "thread id: " + UbSystem::toString(ID));
                  UBLOG(logINFO, "Number of treads: " + UbSystem::toString(omp_get_num_threads()));
               }
               i++;
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
         for (int x3 = 0; x3 < dimensions[2]; x3++)
            for (int x2 = 0; x2 < dimensions[1]; x2++)
               for (int x1 = lstartX1; x1 < lstopX1; x1++)
               {
                  vaFlucVxMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaFlucVyMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaFlucVzMatrix(x1, x2, x3) = receiveBuffer[c++];
                  vaFlucPrMatrix(x1, x2, x3) = receiveBuffer[c++];

                  vaStressXX(x1, x2, x3) = receiveBuffer[c++];
                  vaStressYY(x1, x2, x3) = receiveBuffer[c++];
                  vaStressZZ(x1, x2, x3) = receiveBuffer[c++];
                  vaStressXY(x1, x2, x3) = receiveBuffer[c++];
                  vaStressXZ(x1, x2, x3) = receiveBuffer[c++];
                  vaStressYZ(x1, x2, x3) = receiveBuffer[c++];
               }
      }
   }
   else
   {
      vector<double> sendBuffer;
      for (int x3 = 0; x3 < dimensions[2]; x3++)
         for (int x2 = 0; x2 < dimensions[1]; x2++)
            for (int x1 = startX1; x1 < stopX1; x1++)
            {
               sendBuffer.push_back(vaFlucVxMatrix(x1, x2, x3));
               sendBuffer.push_back(vaFlucVyMatrix(x1, x2, x3));
               sendBuffer.push_back(vaFlucVzMatrix(x1, x2, x3));
               sendBuffer.push_back(vaFlucPrMatrix(x1, x2, x3));

               sendBuffer.push_back(vaStressXX(x1, x2, x3));
               sendBuffer.push_back(vaStressYY(x1, x2, x3));
               sendBuffer.push_back(vaStressZZ(x1, x2, x3));
               sendBuffer.push_back(vaStressXY(x1, x2, x3));
               sendBuffer.push_back(vaStressXZ(x1, x2, x3));
               sendBuffer.push_back(vaStressYZ(x1, x2, x3));
            }
      int count = (int)sendBuffer.size();
      MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&sendBuffer[0], count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&startX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&stopX1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   }

   timer_averaging->StopTimer();
   UBLOG(logINFO, "volume averaging fluct and stress: end");
   UBLOG(logINFO, "volume averaging fluct and stress time: " + UbSystem::toString(timer_averaging->GetElapsedTime()) + " s");
}
void Averaging::writeVolumeAveragingValuesToBinaryFiles(std::string ffname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(vaVxMatrix, ffname + "Vx" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaVyMatrix, ffname + "Vy" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaVzMatrix, ffname + "Vz" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaPrMatrix, ffname + "Pr" + UbSystem::toString(timeStep) + ".bin");
}
void Averaging::readVolumeAveragingValuesFromBinaryFiles(std::string fname, int timeStep)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + UbSystem::toString(timeStep) + ".bin", vaVxMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + UbSystem::toString(timeStep) + ".bin", vaVyMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + UbSystem::toString(timeStep) + ".bin", vaVzMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + UbSystem::toString(timeStep) + ".bin", vaPrMatrix);
}
void Averaging::writeMeanVolumeAveragingValuesToBinaryFiles(std::string ffname)
{
   writeMatrixToBinaryFiles<double>(vaMeanVxMatrix, ffname + "Vx" + ".bin");
   writeMatrixToBinaryFiles<double>(vaMeanVyMatrix, ffname + "Vy" + ".bin");
   writeMatrixToBinaryFiles<double>(vaMeanVzMatrix, ffname + "Vz" + ".bin");
   writeMatrixToBinaryFiles<double>(vaMeanPrMatrix, ffname + "Pr" + ".bin");
}
void Averaging::readMeanVolumeAveragingValuesFromBinaryFiles(std::string ffname)
{
   readMatrixFromBinaryFiles<double>(ffname + "Vx" + ".bin", vaMeanVxMatrix);
   readMatrixFromBinaryFiles<double>(ffname + "Vy" + ".bin", vaMeanVyMatrix);
   readMatrixFromBinaryFiles<double>(ffname + "Vz" + ".bin", vaMeanVzMatrix);
   readMatrixFromBinaryFiles<double>(ffname + "Pr" + ".bin", vaMeanPrMatrix);
}
//void Averaging::readVolumeAveragingFluctStressValuesFromBinaryFiles(std::string fname, int timeStep)
//{
//   readMatrixFromBinaryFiles<double>(fname + "fluctVx" + UbSystem::toString(timeStep) + ".bin", vaFlucVxMatrix);
//   readMatrixFromBinaryFiles<double>(fname + "fluctVy" + UbSystem::toString(timeStep) + ".bin", vaFlucVyMatrix);
//   readMatrixFromBinaryFiles<double>(fname + "fluctVz" + UbSystem::toString(timeStep) + ".bin", vaFlucVzMatrix);
//   readMatrixFromBinaryFiles<double>(fname + "fluctPr" + UbSystem::toString(timeStep) + ".bin", vaFlucPrMatrix);
//
//   readMatrixFromBinaryFiles<double>(fname + "stressXX" + UbSystem::toString(timeStep) + ".bin", vaStressXX);
//   readMatrixFromBinaryFiles<double>(fname + "stressYY" + UbSystem::toString(timeStep) + ".bin", vaStressYY);
//   readMatrixFromBinaryFiles<double>(fname + "stressZZ" + UbSystem::toString(timeStep) + ".bin", vaStressZZ);
//   readMatrixFromBinaryFiles<double>(fname + "stressXY" + UbSystem::toString(timeStep) + ".bin", vaStressXY);
//   readMatrixFromBinaryFiles<double>(fname + "stressXZ" + UbSystem::toString(timeStep) + ".bin", vaStressXZ);
//   readMatrixFromBinaryFiles<double>(fname + "stressYZ" + UbSystem::toString(timeStep) + ".bin", vaStressYZ);
//}

//------------------------------ fluctuations -----------------------
void Averaging::initFluctuations()
{
   flucVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   flucVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   flucVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   flucPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVxMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVyMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucVzMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaFlucPrMatrix.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initSumOfVaFluctuations()
{
   sumVaFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanOfVaFluctuations()
{
   meanVaFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::fluctuationsStress()
{
   vector<double>& vxF = vxMatrix.getDataVector();
   vector<double>& vyF = vyMatrix.getDataVector();
   vector<double>& vzF = vzMatrix.getDataVector();
   vector<double>& prF = prMatrix.getDataVector();

   vector<double>& vxMean = vaMeanVxMatrix.getDataVector();
   vector<double>& vyMean = vaMeanVyMatrix.getDataVector();
   vector<double>& vzMean = vaMeanVzMatrix.getDataVector();
   vector<double>& prMean = vaMeanPrMatrix.getDataVector();

   vector<double>& vxFluc = flucVxMatrix.getDataVector();
   vector<double>& vyFluc = flucVyMatrix.getDataVector();
   vector<double>& vzFluc = flucVzMatrix.getDataVector();
   vector<double>& prFluc = flucPrMatrix.getDataVector();

   vector<double>& XXStress = StressXX.getDataVector();
   vector<double>& YYStress = StressYY.getDataVector();
   vector<double>& ZZStress = StressZZ.getDataVector();
   vector<double>& XYStress = StressXY.getDataVector();
   vector<double>& XZStress = StressXZ.getDataVector();
   vector<double>& YZStress = StressYZ.getDataVector();

   int size = (int)vxF.size();

   for (int i = 0; i < size; i++)
   {
      vxFluc[i] = vxF[i] - vxMean[i];
      vyFluc[i] = vyF[i] - vyMean[i];
      vzFluc[i] = vzF[i] - vzMean[i];
      prFluc[i] = prF[i] - prMean[i];

      XXStress[i] = vxFluc[i] * vxFluc[i];
      YYStress[i] = vyFluc[i] * vyFluc[i];
      ZZStress[i] = vzFluc[i] * vzFluc[i];
      XYStress[i] = vxFluc[i] * vyFluc[i];
      XZStress[i] = vxFluc[i] * vzFluc[i];
      YZStress[i] = vyFluc[i] * vzFluc[i];
   }
}
void Averaging::fluctuationsStress2()
{
   vector<double>& vxVa = vaVxMatrix.getDataVector();
   vector<double>& vyVa = vaVyMatrix.getDataVector();
   vector<double>& vzVa = vaVzMatrix.getDataVector();
   vector<double>& prVa = vaPrMatrix.getDataVector();

   vector<double>& vxMean = vaMeanVxMatrix.getDataVector();
   vector<double>& vyMean = vaMeanVyMatrix.getDataVector();
   vector<double>& vzMean = vaMeanVzMatrix.getDataVector();
   vector<double>& prMean = vaMeanPrMatrix.getDataVector();

   vector<double>& vxFluc = vaFlucVxMatrix.getDataVector();
   vector<double>& vyFluc = vaFlucVyMatrix.getDataVector();
   vector<double>& vzFluc = vaFlucVzMatrix.getDataVector();
   vector<double>& prFluc = vaFlucPrMatrix.getDataVector();

   vector<double>& XXStress = vaStressXX.getDataVector();
   vector<double>& YYStress = vaStressYY.getDataVector();
   vector<double>& ZZStress = vaStressZZ.getDataVector();
   vector<double>& XYStress = vaStressXY.getDataVector();
   vector<double>& XZStress = vaStressXZ.getDataVector();
   vector<double>& YZStress = vaStressYZ.getDataVector();

   int size = (int)vxVa.size();

   for (int i = 0; i < size; i++)
   {
      vxFluc[i] = vxVa[i] - vxMean[i];
      vyFluc[i] = vyVa[i] - vyMean[i];
      vzFluc[i] = vzVa[i] - vzMean[i];
      prFluc[i] = prVa[i] - prMean[i];

      XXStress[i] = vxFluc[i] * vxFluc[i];
      YYStress[i] = vyFluc[i] * vyFluc[i];
      ZZStress[i] = vzFluc[i] * vzFluc[i];
      XYStress[i] = vxFluc[i] * vyFluc[i];
      XZStress[i] = vxFluc[i] * vzFluc[i];
      YZStress[i] = vyFluc[i] * vzFluc[i];
   }
}


void Averaging::sumOfVaFluctuations()
{
   static int counter = 0;
   vector<double>& vxFluc = vaFlucVxMatrix.getDataVector();
   vector<double>& vyFluc = vaFlucVyMatrix.getDataVector();
   vector<double>& vzFluc = vaFlucVzMatrix.getDataVector();
   vector<double>& prFluc = vaFlucPrMatrix.getDataVector();

   vector<double>& SumFlucVx = sumVaFlucVx.getDataVector();
   vector<double>& SumFlucVy = sumVaFlucVy.getDataVector();
   vector<double>& SumFlucVz = sumVaFlucVz.getDataVector();
   vector<double>& SumFlucPr = sumVaFlucPr.getDataVector();

   int size = (int)vxFluc.size();

   for (int i = 0; i < size; i++)
   {
      SumFlucVx[i] += vxFluc[i];
      SumFlucVy[i] += vyFluc[i];
      SumFlucVz[i] += vzFluc[i];
      SumFlucPr[i] += prFluc[i];
      counter++;
   }

}
void Averaging::meanOfVaFluctuations(int numberOfTimeSteps)
{
   vector<double>& MeanFlucVx = meanVaFlucVx.getDataVector();
   vector<double>& MeanFlucVy = meanVaFlucVy.getDataVector();
   vector<double>& MeanFlucVz = meanVaFlucVz.getDataVector();
   vector<double>& MeanFlucPr = meanVaFlucPr.getDataVector();

   vector<double>& SumFlucVx = sumVaFlucVx.getDataVector();
   vector<double>& SumFlucVy = sumVaFlucVy.getDataVector();
   vector<double>& SumFlucVz = sumVaFlucVz.getDataVector();
   vector<double>& SumFlucPr = sumVaFlucPr.getDataVector();

   int size = (int)SumFlucVx.size();

   for (int i = 0; i < size; i++)
   {
      MeanFlucVx[i] = SumFlucVx[i] / numberOfTimeSteps;
      MeanFlucVy[i] = SumFlucVy[i] / numberOfTimeSteps;
      MeanFlucVz[i] = SumFlucVz[i] / numberOfTimeSteps;
      MeanFlucPr[i] = SumFlucPr[i] / numberOfTimeSteps;
   }
}
void Averaging::writeVaFluctuationsToBinaryFiles(std::string fname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(vaFlucVxMatrix, fname + "Vx" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaFlucVyMatrix, fname + "Vy" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaFlucVzMatrix, fname + "Vz" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaFlucPrMatrix, fname + "Pr" + UbSystem::toString(timeStep) + ".bin");
}
void Averaging::readVaFluctuationsFromBinaryFiles(std::string fname, int timeStep)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + UbSystem::toString(timeStep) + ".bin", vaFlucVxMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + UbSystem::toString(timeStep) + ".bin", vaFlucVyMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + UbSystem::toString(timeStep) + ".bin", vaFlucVzMatrix);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + UbSystem::toString(timeStep) + ".bin", vaFlucPrMatrix);
}
void Averaging::initMeanOfVolumeAveragedValues()
{
   meanVaFlucVx.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucVy.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucVz.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaFlucPr.resize(dimensions[0], dimensions[1], dimensions[2], 0);

   meanVaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::sumVolumeAveragedValues()
{
   static int counter = 0;
   vector<double>& vxFluc = vaFlucVxMatrix.getDataVector();
   vector<double>& vyFluc = vaFlucVyMatrix.getDataVector();
   vector<double>& vzFluc = vaFlucVzMatrix.getDataVector();
   vector<double>& prFluc = vaFlucPrMatrix.getDataVector();

   vector<double>& sumFlucVx = meanVaFlucVx.getDataVector();
   vector<double>& sumFlucVy = meanVaFlucVy.getDataVector();
   vector<double>& sumFlucVz = meanVaFlucVz.getDataVector();
   vector<double>& sumFlucPr = meanVaFlucPr.getDataVector();

   vector<double>& sumStressXX = meanVaStressXX.getDataVector();
   vector<double>& sumStressYY = meanVaStressYY.getDataVector();
   vector<double>& sumStressZZ = meanVaStressZZ.getDataVector();
   vector<double>& sumStressXY = meanVaStressXY.getDataVector();
   vector<double>& sumStressXZ = meanVaStressXZ.getDataVector();
   vector<double>& sumStressYZ = meanVaStressYZ.getDataVector();

   int size = (int)vxFluc.size();

   for (int i = 0; i < size; i++)
   {
      sumFlucVx[i] += vxFluc[i];
      sumFlucVy[i] += vyFluc[i];
      sumFlucVz[i] += vzFluc[i];
      sumFlucPr[i] += prFluc[i];

      sumStressXX[i] += vxFluc[i] * vxFluc[i];
      sumStressYY[i] += vyFluc[i] * vyFluc[i];
      sumStressZZ[i] += vzFluc[i] * vzFluc[i];
      sumStressXY[i] += vxFluc[i] * vyFluc[i];
      sumStressXZ[i] += vxFluc[i] * vzFluc[i];
      sumStressYZ[i] += vyFluc[i] * vzFluc[i];

      counter++;
   }
}
void Averaging::computeVolumeAveragedValues(int numberOfTimeSteps)
{
   vector<double>& meanFlucVx = meanVaFlucVx.getDataVector();
   vector<double>& meanFlucVy = meanVaFlucVy.getDataVector();
   vector<double>& meanFlucVz = meanVaFlucVz.getDataVector();
   vector<double>& meanFlucPr = meanVaFlucPr.getDataVector();

   vector<double>& meanStressXX = meanVaStressXX.getDataVector();
   vector<double>& meanStressYY = meanVaStressYY.getDataVector();
   vector<double>& meanStressZZ = meanVaStressZZ.getDataVector();
   vector<double>& meanStressXY = meanVaStressXY.getDataVector();
   vector<double>& meanStressXZ = meanVaStressXZ.getDataVector();
   vector<double>& meanStressYZ = meanVaStressYZ.getDataVector();

   int size = (int)meanFlucVx.size();

   for (int i = 0; i < size; i++)
   {
      meanFlucVx[i] = meanFlucVx[i] / numberOfTimeSteps;
      meanFlucVy[i] = meanFlucVy[i] / numberOfTimeSteps;
      meanFlucVz[i] = meanFlucVz[i] / numberOfTimeSteps;
      meanFlucPr[i] = meanFlucPr[i] / numberOfTimeSteps;

      meanStressXX[i] = meanStressXX[i] / numberOfTimeSteps;
      meanStressYY[i] = meanStressYY[i] / numberOfTimeSteps;
      meanStressZZ[i] = meanStressZZ[i] / numberOfTimeSteps;
      meanStressXY[i] = meanStressXY[i] / numberOfTimeSteps;
      meanStressXZ[i] = meanStressXZ[i] / numberOfTimeSteps;
      meanStressYZ[i] = meanStressYZ[i] / numberOfTimeSteps;
   }
}
void Averaging::writeVolumeAveragedValuesToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<double>(meanVaFlucVx, fname + "Vx" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaFlucVy, fname + "Vy" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaFlucVz, fname + "Vz" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaFlucPr, fname + "Pr" + ".bin");

   writeMatrixToBinaryFiles<double>(meanVaStressXX, fname + "XX" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressYY, fname + "YY" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressZZ, fname + "ZZ" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressXY, fname + "XY" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressXZ, fname + "XZ" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressYZ, fname + "YZ" + ".bin");
}
void Averaging::readVolumeAveragedValuesFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + ".bin", meanVaFlucVx);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + ".bin", meanVaFlucVy);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + ".bin", meanVaFlucVz);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + ".bin", meanVaFlucPr);

   readMatrixFromBinaryFiles<double>(fname + "XX" + ".bin", meanVaStressXX);
   readMatrixFromBinaryFiles<double>(fname + "YY" + ".bin", meanVaStressYY);
   readMatrixFromBinaryFiles<double>(fname + "ZZ" + ".bin", meanVaStressZZ);
   readMatrixFromBinaryFiles<double>(fname + "XY" + ".bin", meanVaStressXY);
   readMatrixFromBinaryFiles<double>(fname + "XZ" + ".bin", meanVaStressXZ);
   readMatrixFromBinaryFiles<double>(fname + "YZ" + ".bin", meanVaStressYZ);
}
void Averaging::writeMeanVaFluctuationsToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<double>(meanVaFlucVx, fname + "Vx" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaFlucVy, fname + "Vy" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaFlucVz, fname + "Vz" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaFlucPr, fname + "Pr" + ".bin");
}
void Averaging::readMeanVaFluctuationsFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<double>(fname + "Vx" + ".bin", meanVaFlucVx);
   readMatrixFromBinaryFiles<double>(fname + "Vy" + ".bin", meanVaFlucVy);
   readMatrixFromBinaryFiles<double>(fname + "Vz" + ".bin", meanVaFlucVz);
   readMatrixFromBinaryFiles<double>(fname + "Pr" + ".bin", meanVaFlucPr);
}
void Averaging::writeMeanOfVaFluctuationsToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { meanVaFlucVx, meanVaFlucVy, meanVaFlucVz, meanVaFlucPr };
   writeMatrixToImageFile(output, matrix);
}
void Averaging::writeFluctuationsToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { flucVxMatrix, flucVyMatrix, flucVzMatrix, flucPrMatrix };
   writeMatrixToImageFile(output, matrix);
}
void Averaging::writeVaFluctuationsToImageFile(std::string output)
{
   array < CbArray3D<double>, 4 > matrix = { vaFlucVxMatrix, vaFlucVyMatrix, vaFlucVzMatrix, vaFlucPrMatrix };
   writeMatrixToImageFile(output, matrix);
}

//----------------------------- stress -----------------------------
void Averaging::initStresses()
{
   StressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   StressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   vaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initSumOfVaStresses()
{
   sumVaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   sumVaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::initMeanOfVaStresses()
{
   meanVaStressXX.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressYY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressZZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressXY.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressXZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
   meanVaStressYZ.resize(dimensions[0], dimensions[1], dimensions[2], 0);
}
void Averaging::sumOfVaStresses()
{
   vector<double>& XXStress = vaStressXX.getDataVector();
   vector<double>& YYStress = vaStressYY.getDataVector();
   vector<double>& ZZStress = vaStressZZ.getDataVector();
   vector<double>& XYStress = vaStressXY.getDataVector();
   vector<double>& XZStress = vaStressXZ.getDataVector();
   vector<double>& YZStress = vaStressYZ.getDataVector();

   vector<double>& XXSum = sumVaStressXX.getDataVector();
   vector<double>& YYSum = sumVaStressYY.getDataVector();
   vector<double>& ZZSum = sumVaStressZZ.getDataVector();
   vector<double>& XYSum = sumVaStressXY.getDataVector();
   vector<double>& XZSum = sumVaStressXZ.getDataVector();
   vector<double>& YZSum = sumVaStressYZ.getDataVector();
                           
   int size = (int)XXStress.size();

   for (int i = 0; i < size; i++)
   {
      XXSum[i] += XXStress[i];
      YYSum[i] += YYStress[i];
      ZZSum[i] += ZZStress[i];
      XYSum[i] += XYStress[i];
      XZSum[i] += XZStress[i];
      YZSum[i] += YZStress[i];
   }
}
void Averaging::meanOfVaStresses(int numberOfTimeSteps)
{
   vector<double>& XXSum = sumVaStressXX.getDataVector();
   vector<double>& YYSum = sumVaStressYY.getDataVector();
   vector<double>& ZZSum = sumVaStressZZ.getDataVector();
   vector<double>& XYSum = sumVaStressXY.getDataVector();
   vector<double>& XZSum = sumVaStressXZ.getDataVector();
   vector<double>& YZSum = sumVaStressYZ.getDataVector();

   vector<double>& XXMean = meanVaStressXX.getDataVector();
   vector<double>& YYMean = meanVaStressYY.getDataVector();
   vector<double>& ZZMean = meanVaStressZZ.getDataVector();
   vector<double>& XYMean = meanVaStressXY.getDataVector();
   vector<double>& XZMean = meanVaStressXZ.getDataVector();
   vector<double>& YZMean = meanVaStressYZ.getDataVector();

   int size = (int)XXSum.size();

   for (int i = 0; i < size; i++)
   {
      XXMean[i] = XXSum[i] / numberOfTimeSteps;
      YYMean[i] = YYSum[i] / numberOfTimeSteps;
      ZZMean[i] = ZZSum[i] / numberOfTimeSteps;
      XYMean[i] = XYSum[i] / numberOfTimeSteps;
      XZMean[i] = XZSum[i] / numberOfTimeSteps;
      YZMean[i] = YZSum[i] / numberOfTimeSteps;
   }
}
void Averaging::writeVaStressesToBinaryFiles(std::string fname, int timeStep)
{
   writeMatrixToBinaryFiles<double>(vaStressXX, fname + "XX" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaStressYY, fname + "YY" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaStressZZ, fname + "ZZ" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaStressXY, fname + "XY" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaStressXZ, fname + "XZ" + UbSystem::toString(timeStep) + ".bin");
   writeMatrixToBinaryFiles<double>(vaStressYZ, fname + "YZ" + UbSystem::toString(timeStep) + ".bin");
}
void Averaging::readVaStressesFromBinaryFiles(std::string fname, int timeStep)
{
   readMatrixFromBinaryFiles<double>(fname + "XX" + UbSystem::toString(timeStep) + ".bin", vaStressXX);
   readMatrixFromBinaryFiles<double>(fname + "YY" + UbSystem::toString(timeStep) + ".bin", vaStressYY);
   readMatrixFromBinaryFiles<double>(fname + "ZZ" + UbSystem::toString(timeStep) + ".bin", vaStressZZ);
   readMatrixFromBinaryFiles<double>(fname + "XY" + UbSystem::toString(timeStep) + ".bin", vaStressXY);
   readMatrixFromBinaryFiles<double>(fname + "XZ" + UbSystem::toString(timeStep) + ".bin", vaStressXZ);
   readMatrixFromBinaryFiles<double>(fname + "YZ" + UbSystem::toString(timeStep) + ".bin", vaStressYZ);
}
void Averaging::writeMeanVaStressesToBinaryFiles(std::string fname)
{
   writeMatrixToBinaryFiles<double>(meanVaStressXX, fname + "XX" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressYY, fname + "YY" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressZZ, fname + "ZZ" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressXY, fname + "XY" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressXZ, fname + "XZ" + ".bin");
   writeMatrixToBinaryFiles<double>(meanVaStressYZ, fname + "YZ" + ".bin");
}
void Averaging::readMeanVaStressesFromBinaryFiles(std::string fname)
{
   readMatrixFromBinaryFiles<double>(fname + "XX" + ".bin", meanVaStressXX);
   readMatrixFromBinaryFiles<double>(fname + "YY" + ".bin", meanVaStressYY);
   readMatrixFromBinaryFiles<double>(fname + "ZZ" + ".bin", meanVaStressZZ);
   readMatrixFromBinaryFiles<double>(fname + "XY" + ".bin", meanVaStressXY);
   readMatrixFromBinaryFiles<double>(fname + "XZ" + ".bin", meanVaStressXZ);
   readMatrixFromBinaryFiles<double>(fname + "YZ" + ".bin", meanVaStressYZ);
}

void Averaging::writeMeanOfVaStressesToImageFile(std::string ffname)
{
   array < CbArray3D<double>, 4 > matrix = { meanVaStressXX, meanVaStressYY, meanVaStressZZ, meanVaStressXY };
   writeMatrixToImageFile(ffname, matrix);
}

//------------------------------------ planar --------------------------
void Averaging::initPlanarAveraging()
{
   PlanarVx.resize(dimensions[2], 0);
   PlanarVy.resize(dimensions[2], 0);
   PlanarVz.resize(dimensions[2], 0);
   PlanarPr.resize(dimensions[2], 0);

   PlanarFlucVx.resize(dimensions[2], 0);
   PlanarFlucVy.resize(dimensions[2], 0);
   PlanarFlucVz.resize(dimensions[2], 0);
   PlanarFlucPr.resize(dimensions[2], 0);

   PlanarStressXX.resize(dimensions[2], 0);
   PlanarStressYY.resize(dimensions[2], 0);
   PlanarStressZZ.resize(dimensions[2], 0);
   PlanarStressXY.resize(dimensions[2], 0);
   PlanarStressXZ.resize(dimensions[2], 0);
   PlanarStressYZ.resize(dimensions[2], 0);
}
void Averaging::planarAveraging()
{
   double numberof_XY_points = (double)dimensions[0] * (double)dimensions[1];

   for (int z = 0; z < dimensions[2]; z++)
   {
      double sumVx = 0, sumVy = 0, sumVz = 0, sumPr = 0;
      double sumFluctVx = 0, sumFluctVy = 0, sumFluctVz = 0, sumFluctPr = 0;
      double sumStressXX = 0, sumStressYY = 0, sumStressZZ = 0, sumStressXY = 0, sumStressXZ = 0, sumStressYZ = 0;
      for (int y = 0; y < dimensions[1]; y++)
         for (int x = 0; x < dimensions[0]; x++)
         {
            sumVx += vaMeanVxMatrix(x, y, z);
            sumVy += vaMeanVyMatrix(x, y, z);
            sumVz += vaMeanVzMatrix(x, y, z);
            sumPr += vaMeanPrMatrix(x, y, z);

            sumFluctVx += meanVaFlucVx(x, y, z);
            sumFluctVy += meanVaFlucVy(x, y, z);
            sumFluctVz += meanVaFlucVz(x, y, z);
            sumFluctPr += meanVaFlucPr(x, y, z);

            sumStressXX += meanVaStressXX(x, y, z);
            sumStressYY += meanVaStressYY(x, y, z);
            sumStressZZ += meanVaStressZZ(x, y, z);
            sumStressXY += meanVaStressXY(x, y, z);
            sumStressXZ += meanVaStressXZ(x, y, z);
            sumStressYZ += meanVaStressYZ(x, y, z);
         }
      PlanarVx[z] = sumVx / numberof_XY_points;
      PlanarVy[z] = sumVy / numberof_XY_points;
      PlanarVz[z] = sumVz / numberof_XY_points;
      PlanarPr[z] = sumPr / numberof_XY_points;

      PlanarFlucVx[z] = sumFluctVx / numberof_XY_points;
      PlanarFlucVy[z] = sumFluctVy / numberof_XY_points;
      PlanarFlucVz[z] = sumFluctVz / numberof_XY_points;
      PlanarFlucPr[z] = sumFluctPr / numberof_XY_points;

      PlanarStressXX[z] = sumStressXX / numberof_XY_points;
      PlanarStressYY[z] = sumStressYY / numberof_XY_points;
      PlanarStressZZ[z] = sumStressZZ / numberof_XY_points;
      PlanarStressXY[z] = sumStressXY / numberof_XY_points;
      PlanarStressXZ[z] = sumStressXZ / numberof_XY_points;
      PlanarStressYZ[z] = sumStressYZ / numberof_XY_points;
   }
}
   
void Averaging::writeToCSV(std::string path, double origin, double deltax)
   {
      std::ofstream ostr;
      std::string fname = path + ".csv";

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
      ostr << "z;Vx;Vy;Vz;Pr;FlucVx;FlucVy;FlucVz;FlucPr;StressXX;StressYY;StressZZ;StressXY;StressXZ;StressYZ\n";
      for (int i = 0; i < dimensions[2]; i++)
      {
         double z = origin + (deltax*i);
         ostr << z << ";" << PlanarVx[i] << ";" << PlanarVy[i] << ";" << PlanarVz[i] << ";" << PlanarPr[i] << ";"  << PlanarFlucVx[i] << ";" << PlanarFlucVy[i] << ";" << PlanarFlucVz[i] << ";" << PlanarFlucPr[i] << ";" << PlanarStressXX[i] << ";" << PlanarStressYY[i] << ";" << PlanarStressZZ[i] << ";" << PlanarStressXY[i] << ";" << PlanarStressXZ[i] << ";" << PlanarStressYZ[i] << "\n";
      }
      ostr.close();
   }

   void Averaging::writeToCSV2(std::string path, double origin, double deltax)
   {
      std::ofstream ostr;
      std::string fname = path + ".csv";

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
      ostr << "z;Vx;Vy;Vz;Pr;StressXX;StressYY;StressZZ;StressXY;StressXZ;StressYZ\n";
      for (int i = 0; i < dimensions[2]; i++)
      {
         double z = origin + (deltax*i);
         ostr << z << ";" << PlanarVx[i] << ";" << PlanarVy[i] << ";" << PlanarVz[i] << ";" << PlanarPr[i] << ";"<< PlanarStressXX[i] << ";" << PlanarStressYY[i] << ";" << PlanarStressZZ[i] << ";" << PlanarStressXY[i] << ";" << PlanarStressXZ[i] << ";" << PlanarStressYZ[i] << "\n";
      }
      ostr.close();
   }