#include <iostream>
#include <string>
#include "VirtualFluids.h"

#include "Averaging.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   std::array<int,3> dimensions = {60,40,40};
   Averaging av;
   av.setDimensions(dimensions);
   int geo_extent[6] = { 0, dimensions[0]-1, 0, dimensions[1]-1, 0,dimensions[2]-1 };
   double geo_origin[3] = { 0.5,0.5,0.5 };
   double geo_spacing[3] = { 1,1,1 };
   double deltax = 1;
   av.createGeoMatrix("e:/temp/BreugemChannelAnisotrop/bc/bc0.pvtu", deltax, geo_origin);
   av.writeGeoMatrixToBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/geomatrix.bin");
   av.readGeoMatrixFromBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/geomatrix.bin");
   av.writeGeoMatrixToImageFile("e:/temp/BreugemChannelAnisotrop/va/geoMatrix", geo_extent, geo_origin, geo_spacing);

   av.readGeoMatrix("e:/temp/BreugemChannelAnisotrop/va/geoMatrix.vti");
   av.createMQMatrix("e:/temp/BreugemChannelAnisotrop/mq/mq100.pvtu", deltax, geo_origin);
   av.writeMQMatrixToImageFile("e:/temp/BreugemChannelAnisotrop/va/mqMatrix", geo_extent, geo_origin, geo_spacing);
   return 0;
}
