#include <array>
#include "CbArray3D.h"
#include "UbSystem.h"



class Averaging
{
public:
   void createGeoMatrix(std::string dataNameG, double deltax, double geo_origin[3]);
   void writeGeoMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3]);
   void createMQMatrix(std::string dataNameMQ, double deltax, double geo_origin[3]);
   void writeMQMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3]);
   void averagingWithMPI(double l_real, double l);
   void Averaging::readGeoMatrix(std::string dataNameG);
   void writeGeoMatrixToBinaryFiles(std::string fname);
   void readGeoMatrixFromBinaryFiles(std::string fname);


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

   CbArray3D<double> vaVxMatrix;
   CbArray3D<double> vaVyMatrix;
   CbArray3D<double> vaVzMatrix;
   CbArray3D<double> vaPrMatrix;
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