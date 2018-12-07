#ifndef Postprocessing_h
#define Postprocessing_h

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkTimerLog.h>

namespace VAP
{
   void volumeAveragingWithProbe();
}

namespace VAI
{
   void volumeAveragingWithInterpolation(int n);
}

namespace AA
{
   void AvAllRun(char c);
}

namespace ATAV
{
   void AvAllRun(char c);
}


namespace PP
{
   void CalculateFluctuations(std::vector <std::string> dataNames, std::string meanVelocityData, std::string output);
   void CalculateMeanVelocity(std::vector <std::string> dataName, std::string output);
}


void MarchingCubes();

void StlToVtu(std::string stlfile, std::string vtifile, double spacing[3]);

void VolumeAveragingWithVector();

vtkDataSet* ReadDataSet(std::string fileName);
void SurfaceAveragingMeanVelocity(std::string dataName, double step, std::string output);
void CalculateBulkVelocity(std::string dataName, double origin[3]);
void VolumeAveraging(std::string dataName, double radius, std::string output);

double G(double x, double x_min, double x_max);
double m(double x, double y, double z, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double norm);
void VolumenFilter();
const std::string currentDateTime();
void print(const std::string& s);
void VolumeAveragingBreugem(std::string dataName, double L, double dx, std::string output);
void VolumeAveragingToImage(std::string dataName, double L, double dx, double image_dx, std::string output);
void VolumeAveragingWithMatrix(std::string dataName, std::string dataNameL, std::string dataNameRMS, int numberOfLevels, /*int msize[3],*/ double L, double deltax, std::string output);
void CalculateRMS(std::vector <std::string> dataNames, std::string meanVelocityData, std::string output);
void SurfaceAveragingRMS(std::string dataName, double step, std::string output);
int cint(double x);

//////////////////////////////////////////////////////////////////////////
template<class TReader> vtkDataSet *ReadAnXMLFile(std::string fileName)
{
   vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
   reader->SetFileName(fileName.c_str());
   reader->Update();
   reader->GetOutput()->Register(reader);
   return vtkDataSet::SafeDownCast(reader->GetOutput());
}
//////////////////////////////////////////////////////////////////////////
template<class T> std::string toString(const T& t)
{
   std::ostringstream stream;
   stream << t;
   return stream.str();
}
//////////////////////////////////////////////////////////////////////////
template<class T> T fromString(const std::string& s)
{
   //boolean hack
   if(s == "true")
      return true;
   else if(s == "false")
      return false;
   //////////////
   std::istringstream stream (s);
   T t;
   stream >> t;
   return t;
}
//////////////////////////////////////////////////////////////////////////
template<class T> void writeMatrixToBinaryFiles(std::vector<T> &matrix, std::string fname)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   print("write matrix to " + fname + ": start");
   timer_write->StartTimer();


   FILE *file;
   file = fopen(fname.c_str(), "wb");

   if (file == NULL)
   {
      print("can not open " + fname);
      return;
   }

   fwrite(&matrix[0], sizeof(T), matrix.size(), file);

   fclose(file);

   print("write matrix: end");
   timer_write->StopTimer();
   print("write matrix time: " + toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
template<class T> void readMatrixFromBinaryFiles(std::string fname, std::vector<T> &matrix)
{
   vtkSmartPointer<vtkTimerLog> timer_write = vtkSmartPointer<vtkTimerLog>::New();

   print("read matrix from " + fname + ": start");
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

   print("read matrix: end");
   timer_write->StopTimer();
   print("read matrix time: " + toString(timer_write->GetElapsedTime()) + " s");
}
//////////////////////////////////////////////////////////////////////////
#endif

