#include <array>
#include "CbArray3D.h"




class Averaging
{
public:
   void createGeoMatrix(std::string dataNameG, double deltax, double geo_origin[3]);
   void writeGeoMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3]);
   void createMQMatrix(std::string dataNameMQ, double deltax, double geo_origin[3]);
   void writeMQMatrixToImageFile(std::string output, int geo_extent[6], double geo_origin[3], double geo_spacing[3]);

   std::array<int, 3> getDimensions() const { return dimensions; }
   void setDimensions(std::array<int, 3> val) { dimensions = val; }
protected:
   void getNodeIndexes(double x[3], int ix[3], double origin[3], double deltax);
   double G(double x);
private:
   std::array<int,3> dimensions;
   CbArray3D<int> geoMatrix;

   CbArray3D<double> vxMatrix;
   CbArray3D<double> vyMatrix;
   CbArray3D<double> vzMatrix;
   CbArray3D<double> prMatrix;
};
