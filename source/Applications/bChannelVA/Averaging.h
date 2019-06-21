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
   void readGeoMatrix(std::string dataNameG);
   void writeGeoMatrixToBinaryFiles(std::string fname);
   void readGeoMatrixFromBinaryFiles(std::string fname);

   void createMQMatrix(std::string dataNameMQ);
   void writeMqMatrixToBinaryFiles(std::string fname, int timeStep);
   void readMqMatrixFromBinaryFiles(std::string fname, int timeStep);
   void writeMqMatrixToImageFile(std::string output);
   void writeVaMatrixToImageFile(std::string output);
   void writeVaSumMatrixToImageFile(std::string output);
   void writeMeanMatrixToImageFile(std::string output);
   void writeMatrixToImageFile(std::string output, std::array<CbArray3D<double>, 4> matrix);

   void initVolumeAveragingValues();
   void initVolumeAveragingFluctStressValues();
   void initMeanVolumeAveragingValues();
   void initMeanVolumeAveragingFluctStressValues();
   void volumeAveragingWithMPI(double l_real);
   void volumeAveragingFluctStressWithMPI(double l_real);
   void meanOfVolumeAveragingValues(int numberOfTimeSteps);
   void sumOfVolumeAveragingValues();
   void writeVolumeAveragingValuesToBinaryFiles(std::string ffname, int timeStep);
   void readVolumeAveragingValuesFromBinaryFiles(std::string fname, int timeStep);
   void writeMeanVolumeAveragingValuesToBinaryFiles(std::string ffname);
   void readMeanVolumeAveragingValuesFromBinaryFiles(std::string fname);

   void initFluctuations();
   void initMeanOfVaFluctuations();
   void initSumOfVaFluctuations();
   void fluctuationsStress();
   void meanOfVaFluctuations(int numberOfTimeSteps);
   void sumOfVaFluctuations();
   void writeVaFluctuationsToBinaryFiles(std::string fname, int timeStep);
   void readVaFluctuationsFromBinaryFiles(std::string fname, int timeStep);
   void writeMeanVaFluctuationsToBinaryFiles(std::string ffname);
   void readMeanVaFluctuationsFromBinaryFiles(std::string ffname);
   void writeMeanOfVaFluctuationsToImageFile(std::string ffname);
   void writeFluctuationsToImageFile(std::string ffname);
   void writeVaFluctuationsToImageFile(std::string ffname);

   void initStresses();
   void initSumOfVaStresses();
   void initMeanOfVaStresses();
   void sumOfVaStresses();
   void meanOfVaStresses(int numberOfTimeSteps);
   void writeVaStressesToBinaryFiles(std::string fname, int timeStep);
   void readVaStressesFromBinaryFiles(std::string fname, int timeStep);
   void writeMeanVaStressesToBinaryFiles(std::string ffname);
   void readMeanVaStressesFromBinaryFiles(std::string ffname);

   void initPlanarAveraging();
   void planarAveraging();
 
   void writeToCSV(std::string path, double origin, double deltax);

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

   CbArray3D<double> vaVxMatrix;
   CbArray3D<double> vaVyMatrix;
   CbArray3D<double> vaVzMatrix;
   CbArray3D<double> vaPrMatrix;

   CbArray3D<double> sumVaVxMatrix;
   CbArray3D<double> sumVaVyMatrix;
   CbArray3D<double> sumVaVzMatrix;
   CbArray3D<double> sumVaPrMatrix;

   CbArray3D<double> meanVaVxMatrix;
   CbArray3D<double> meanVaVyMatrix;
   CbArray3D<double> meanVaVzMatrix;
   CbArray3D<double> meanVaPrMatrix;
//----------------------------------------
   CbArray3D<double> FlucVxMatrix;
   CbArray3D<double> FlucVyMatrix;
   CbArray3D<double> FlucVzMatrix;
   CbArray3D<double> FlucPrMatrix;

   CbArray3D<double> vaFlucVxMatrix;
   CbArray3D<double> vaFlucVyMatrix;
   CbArray3D<double> vaFlucVzMatrix;
   CbArray3D<double> vaFlucPrMatrix;

   CbArray3D<double> sumVaFlucVx;
   CbArray3D<double> sumVaFlucVy;
   CbArray3D<double> sumVaFlucVz;
   CbArray3D<double> sumVaFlucPr;

   CbArray3D<double> meanVaFlucVx;
   CbArray3D<double> meanVaFlucVy;
   CbArray3D<double> meanVaFlucVz;
   CbArray3D<double> meanVaFlucPr;
//----------------------------------------
   CbArray3D<double> StressXX;
   CbArray3D<double> StressYY;
   CbArray3D<double> StressZZ;
   CbArray3D<double> StressXY;
   CbArray3D<double> StressXZ;
   CbArray3D<double> StressYZ;

   CbArray3D<double> vaStressXX;
   CbArray3D<double> vaStressYY;
   CbArray3D<double> vaStressZZ;
   CbArray3D<double> vaStressXY;
   CbArray3D<double> vaStressXZ;
   CbArray3D<double> vaStressYZ;

   CbArray3D<double> SumVaStressXX;
   CbArray3D<double> SumVaStressYY;
   CbArray3D<double> SumVaStressZZ;
   CbArray3D<double> SumVaStressXY;
   CbArray3D<double> SumVaStressXZ;
   CbArray3D<double> SumVaStressYZ;

   CbArray3D<double> meanVaStressXX;
   CbArray3D<double> meanVaStressYY;
   CbArray3D<double> meanVaStressZZ;
   CbArray3D<double> meanVaStressXY;
   CbArray3D<double> meanVaStressXZ;
   CbArray3D<double> meanVaStressYZ;
//----------------------------------------
   std::vector<double> PlanarVx;
   std::vector<double> PlanarVy;
   std::vector<double> PlanarVz;
   std::vector<double> PlanarPr;

   std::vector<double> PlanarFlucVx;
   std::vector<double> PlanarFlucVy;
   std::vector<double> PlanarFlucVz;
   std::vector<double> PlanarFlucPr;

   std::vector<double> PlanarStressXX;
   std::vector<double> PlanarStressYY;
   std::vector<double> PlanarStressZZ;
   std::vector<double> PlanarStressXY;
   std::vector<double> PlanarStressXZ;
   std::vector<double> PlanarStressYZ;
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
inline double Averaging::G(double x, double l)
{
   if (fabs(x) <= l)
      return l - fabs(x);
   else
      return 0.0;
}
//////////////////////////////////////////////////////////////////////////