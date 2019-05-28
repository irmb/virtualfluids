#include <array>
#include <fstream>
#include "CbArray3D.h"
#include "UbSystem.h"
#include <vtkTimerLog.h>
#include <vtkSmartPointer.h>


class Averaging
{
public:
   void createGeoMatrix(std::string dataNameG);
   void writeGeoMatrixToImageFile(std::string output);
   void createMQMatrix(std::string dataNameMQ);
   void writeMatrixToImageFile(std::string output, std::array<CbArray3D<double>, 4> matrix);
   void writeMqMatrixToImageFile(std::string output);
   void writeVaMatrixToImageFile(std::string output);
   void writeVaSumMatrixToImageFile(std::string output);
   void writeMeanMatrixToImageFile(std::string output);
   void volumeAveragingWithMPI(double l_real);
   void readGeoMatrix(std::string dataNameG);
   void writeGeoMatrixToBinaryFiles(std::string fname);
   void readGeoMatrixFromBinaryFiles(std::string fname);
   void writeMqMatrixToBinaryFiles(std::string fname, int timeStep);
   void readMqMatrixFromBinaryFiles(std::string fname, int timeStep);
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
   void sumOfStresses();
   void meanOfStresses(int numberOfTimeSteps);
   void planarAveragingMQ(std::array<int, 3> dimensions);
   void writeToCSV(std::string path, double origin, double deltax);
   void readVolumeAveragingValuesFromBinaryFiles(std::string fname, int timeStep);

   std::array<int, 3> getDimensions() const { return dimensions; }
   void setDimensions(std::array<int, 3> val) { dimensions = val; }
   void setExtent(std::array<int, 6> val) { geo_extent = val; }
   void setOrigin(std::array<double, 3> val) { geo_origin = val; }
   void setSpacing(std::array<double, 3> val) { geo_spacing = val; }
   void setDeltaX(double val) { deltax = val; }
protected:
   void getNodeIndexes(std::array<double, 3> x, std::array<int, 3>& ix);
   double G(double x, double l);
   
   template <class T>
   void writeMatrixToBinaryFiles(CbArray3D<T>& matrix, std::string fname);
   template <class T>
   void readMatrixFromBinaryFiles(std::string fname, CbArray3D<T>& matrix);
private:
   std::array<int, 3> dimensions;
   std::array<int, 6> geo_extent;
   std::array<double, 3> geo_origin;
   std::array<double, 3> geo_spacing;
   double deltax;
 
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
template<class T> void Averaging::writeMatrixToBinaryFiles(CbArray3D<T>& matrix, std::string fname)
 {
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   UBLOG(logINFO,"write matrix to " + fname + ": start");
   timer_write->StartTimer();

   std::ofstream ostr;
   ostr.open(fname.c_str(), std::fstream::out | std::fstream::binary);
   
   if (!ostr)
   {
      ostr.clear();
      std::string path = UbSystem::getPathFromString(fname);
      if (path.size() > 0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::fstream::binary); }
      if (!ostr) throw UbException(UB_EXARGS, "couldn't open file " + fname);
   }

   std::vector<T>& vec = matrix.getDataVector();

   ostr.write((char*)& vec[0], sizeof(T)*vec.size());
   ostr.close();

   UBLOG(logINFO,"write matrix: end");
   timer_write->StopTimer();
   UBLOG(logINFO,"write matrix time: " + UbSystem::toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
template<class T> void Averaging::readMatrixFromBinaryFiles(std::string fname, CbArray3D<T>& matrix)
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
   //matrix.resize(lSize);
   matrix.resize(dimensions[0], dimensions[1], dimensions[2]);
   std::vector<T>& vec = matrix.getDataVector();

   if (vec.size() == 0) { fputs("Memory error", stderr); exit(2); }

   // copy the file into the buffer:
   size_t result = fread(&vec[0], sizeof(T), lSize, file);
   if (result != lSize) { fputs("Reading error", stderr); exit(3); }

   fclose(file);

   UBLOG(logINFO,"read matrix: end");
   timer_write->StopTimer();
   UBLOG(logINFO,"read matrix time: " + UbSystem::toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////