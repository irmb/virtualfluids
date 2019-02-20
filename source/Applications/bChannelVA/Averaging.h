#include <array>
#include "CbArray3D.h"
#include "UbSystem.h"
#include <vtkTimerLog.h>
#include <vtkSmartPointer.h>


class Averaging
{
public:
   void createGeoMatrix(std::string dataNameG, double deltax, double geo_origin[3]);
   void writeGeoMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3]);
   void createMQMatrix(std::string dataNameMQ, double deltax, double geo_origin[3]);
   void writeMQMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3]);
   void volumeAveragingWithMPI(double l_real, double l);
   void readGeoMatrix(std::string dataNameG);
   void writeGeoMatrixToBinaryFiles(std::string fname);
   void readGeoMatrixFromBinaryFiles(std::string fname);
   void initVolumeAveragingValues();
   void initMeanVolumeAveragingValues();
   void initFluctuationsofVolumeAveragingValues();
   void initMeanOfFluctuations();
   void initStresses();
   void initSumOfStresses();
   void initMeanOfStresses();
   void initPlanarAveragingMQ();
   void sumOfVolumeAveragingValues();
   void writeVolumeAveragingValuesToBinaryFiles(std::string ffname, int timeStep);
   void meanOfVolumeAveragingValues(int numberOfTimeSteps);
   void writeMeanVolumeAveragingValuesToBinaryFiles(std::string ffname);
   void fluctuationsOfVolumeAveragingValue();
   void sumOfFluctuations();
   void initSumOfFluctuations();
   void writeFluctuationsToBinaryFiles(std::string fname, int timeStep);
   void writeStressesToBinaryFiles(std::string fname, int timeStep);
   void meanOfFluctuations(int numberOfTimeSteps);
   void SumOfStresses();
   void MeanOfStresses(int numberOfTimeSteps);
   void PlanarAveragingMQ(std::array<int, 3> dimensions);
   void WriteToCSV(std::string path, double origin, double deltax);
   void readVolumeAveragingValuesFromBinaryFiles(std::string fname, int timeStep);

   std::array<int, 3> getDimensions() const { return dimensions; }
   void setDimensions(std::array<int, 3> val) { dimensions = val; }
protected:
   void getNodeIndexes(double x[3], int ix[3], double origin[3], double deltax);
   double G(double x, double l);
   
   template <class T>
   void writeMatrixToBinaryFiles(std::vector<T> &matrix, std::string fname);
   template <class T>
   void readMatrixFromBinaryFiles(std::string fname, std::vector<T> &matrix);
private:
   std::array<int,3> dimensions;
   CbArray3D<int> geoMatrix;

   CbArray3D<double> vxMatrix;
   CbArray3D<double> vyMatrix;
   CbArray3D<double> vzMatrix;
   CbArray3D<double> prMatrix;

   CbArray3D<double> sumVaVxMatrix;
   CbArray3D<double> sumVaVyMatrix;
   CbArray3D<double> sumVaVzMatrix;
   CbArray3D<double> sumVaPrMatrix;

   CbArray3D<double> vaVxMatrix;
   CbArray3D<double> vaVyMatrix;
   CbArray3D<double> vaVzMatrix;
   CbArray3D<double> vaPrMatrix;

   CbArray3D<double> meanVxMatrix;
   CbArray3D<double> meanVyMatrix;
   CbArray3D<double> meanVzMatrix;
   CbArray3D<double> meanPrMatrix;

   CbArray3D<double> FlucVxMatrix;
   CbArray3D<double> FlucVyMatrix;
   CbArray3D<double> FlucVzMatrix;
   CbArray3D<double> FlucPrMatrix;

   CbArray3D<double> StressXX;
   CbArray3D<double> StressXY;
   CbArray3D<double> StressXZ;
   CbArray3D<double> StressYX;
   CbArray3D<double> StressYY;
   CbArray3D<double> StressYZ;
   CbArray3D<double> StressZX;
   CbArray3D<double> StressZY;
   CbArray3D<double> StressZZ;

   CbArray3D<double> sumFlucVx;
   CbArray3D<double> sumFlucVy;
   CbArray3D<double> sumFlucVz;
   CbArray3D<double> sumFlucPr;

   CbArray3D<double> meanFlucVx;
   CbArray3D<double> meanFlucVy;
   CbArray3D<double> meanFlucVz;
   CbArray3D<double> meanFlucPr;

   CbArray3D<double> SumStressXX;
   CbArray3D<double> SumStressXY;
   CbArray3D<double> SumStressXZ;
   CbArray3D<double> SumStressYX;
   CbArray3D<double> SumStressYY;
   CbArray3D<double> SumStressYZ;
   CbArray3D<double> SumStressZX;
   CbArray3D<double> SumStressZY;
   CbArray3D<double> SumStressZZ;

   CbArray3D<double> meanStressXX;
   CbArray3D<double> meanStressXY;
   CbArray3D<double> meanStressXZ;
   CbArray3D<double> meanStressYX;
   CbArray3D<double> meanStressYY;
   CbArray3D<double> meanStressYZ;
   CbArray3D<double> meanStressZX;
   CbArray3D<double> meanStressZY;
   CbArray3D<double> meanStressZZ;

   std::vector<double> PlanarVx;
   std::vector<double> PlanarVy;
   std::vector<double> PlanarVz;
   std::vector<double> PlanarPr;

   std::vector<double> PlanarStressXX;
   std::vector<double> PlanarStressXY;
   std::vector<double> PlanarStressXZ;
  
   std::vector<double> PlanarStressYX;
   std::vector<double> PlanarStressYY;
   std::vector<double> PlanarStressYZ;
   
   std::vector<double> PlanarStressZX;
   std::vector<double> PlanarStressZY;
   std::vector<double> PlanarStressZZ;
};

//////////////////////////////////////////////////////////////////////////
template<class T> void Averaging::writeMatrixToBinaryFiles(std::vector<T> &matrix, std::string fname)
 {
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO,"write matrix to " + fname + ": start");
   timer_write->StartTimer();


   FILE *file;
   file = fopen(fname.c_str(), "wb");

   if (file == NULL)
   {
      UBLOG(logINFO,"can not open " + fname);
      return;
   }

   fwrite(&matrix[0], sizeof(T), matrix.size(), file);

   fclose(file);

   UBLOG(logINFO,"write matrix: end");
   timer_write->StopTimer();
   UBLOG(logINFO,"write matrix time: " + UbSystem::toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
template<class T> void Averaging::readMatrixFromBinaryFiles(std::string fname, std::vector<T> &matrix)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO,"read matrix from " + fname + ": start");
   timer_write->StartTimer();

   FILE *file;
   file = fopen(fname.c_str(), "rb");

   if (file==NULL) { fputs("File error", stderr); exit(1); }

   // obtain file size:
   fseek(file, 0, SEEK_END);
   long lSize = ftell(file)/sizeof(T);
   rewind(file);

   // allocate memory to contain the whole file:
   matrix.resize(lSize);
   if (matrix.size() == 0) { fputs("Memory error", stderr); exit(2); }

   // copy the file into the buffer:
   size_t result = fread(&matrix[0], sizeof(T), lSize, file);
   if (result != lSize) { fputs("Reading error", stderr); exit(3); }

   fclose(file);

   UBLOG(logINFO,"read matrix: end");
   timer_write->StopTimer();
   UBLOG(logINFO,"read matrix time: " + UbSystem::toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////