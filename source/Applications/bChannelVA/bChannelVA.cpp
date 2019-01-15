#include <iostream>
#include <string>
#include "VirtualFluids.h"

#include "Averaging.h"

//using namespace std;


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();

   std::array<int, 3> dimensions = { 60,40,40 };
   Averaging av;
   av.setDimensions(dimensions);
   int geo_extent[6] = { 0, dimensions[0] - 1, 0, dimensions[1] - 1, 0,dimensions[2] - 1 };
   double geo_origin[3] = { 0.5,0.5,0.5 };
   double geo_spacing[3] = { 1,1,1 };
   double deltax = 1;
   double real_l = 40;
   double l = 40;
   av.createGeoMatrix("e:/temp/BreugemChannelAnisotrop/bc/bc0.pvtu", deltax, geo_origin);
   av.writeGeoMatrixToBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/geomatrix.bin");
   int numberOfTimeSteps = 5;
   int numberOfGridPoints = dimensions[0]* dimensions[1]* dimensions[2];
   av.initVolumeAveragingValues();
   for (int t = 0; t < numberOfTimeSteps; t++)
   {
      av.createMQMatrix("e:/temp/BreugemChannelAnisotrop/mq/mq"+UbSystem::toString((t+1)*100)+".pvtu", deltax, geo_origin);
      av.volumeAveragingWithMPI(real_l, l);
      av.sumOfVolumeAveragingValues();
      av.writeVolumeAveragingValuesToBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/vav/vav",t);
   }
   av.initMeanVolumeAveragingValues();
   av.meanOfVolumeAveragingValues(numberOfTimeSteps);
   av.writeMeanVolumeAveragingValuesToBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/mean/mean");
   av.initFluctuationsofVolumeAveragingValues();
   av.initSumOfFluctuations();
   av.initMeanOfFluctuations();
   av.initStresses();
   av.initPlanarAveragingMQ();
   for (int i = 0; i < numberOfTimeSteps; i++)
   {
      av.readVolumeAveragingValuesFromBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/vav/vav",i);
      av.fluctuationsOfVolumeAveragingValue();
      av.writeFluctuationsToBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/av/Fluc", i);
      av.writeStressesToBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/av/Stresses", i);
      av.sumOfFluctuations();
      av.SumOfStresses();
   }
   av.meanOfFluctuations(numberOfTimeSteps);
   av.MeanOfStresses(numberOfTimeSteps);
   av.PlanarAveragingMQ(dimensions);
   av.WriteToCSV("e:/temp/BreugemChannelAnisotrop/va/av", geo_origin[2], deltax);
   //      //av.readGeoMatrixFromBinaryFiles("e:/temp/BreugemChannelAnisotrop/va/geomatrix.bin");
   //      //av.writeGeoMatrixToImageFile("e:/temp/BreugemChannelAnisotrop/va/geoMatrix", geo_extent, geo_origin, geo_spacing);

   //      //av.readGeoMatrix("e:/temp/BreugemChannelAnisotrop/va/geoMatrix.vti");
   //      //av.createMQMatrix("e:/temp/BreugemChannelAnisotrop/mq/mq100.pvtu", deltax, geo_origin);
   //      //av.writeMQMatrixToImageFile("e:/temp/BreugemChannelAnisotrop/va/mqMatrix", geo_extent, geo_origin, geo_spacing);
   //      return 0;
   //   }
}