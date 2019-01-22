#include <iostream>
#include <string>
#include "VirtualFluids.h"

#include "Averaging.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   SPtr<Communicator> comm = MPICommunicator::getInstance();

   string pathIn  = "/work/koskuche/BreugemChannelIsotropic";
   string pathOut = "/work/i5042202/BreugemChannelIsotropic";

   std::array<int, 3> dimensions ={ 600,400,400 };
   Averaging av;
   av.setDimensions(dimensions);
   int geo_extent[6] ={ 0, dimensions[0] - 1, 0, dimensions[1] - 1, 0,dimensions[2] - 1 };
   double geo_origin[3] ={ 0.5,0.5,0.5 };
   double geo_spacing[3] ={ 1,1,1 };
   double deltax = 1;
   double real_l = 40;
   double l = 40;
   av.createGeoMatrix(pathIn + "/bc/bc0.pvtu", deltax, geo_origin);
   av.writeGeoMatrixToBinaryFiles(pathOut + "/va/geomatrix.bin");
   int startTimeStep = 600000;
   int timeStep = 10000;
   int numberOfTimeSteps = 600000;
   int numberOfGridPoints = dimensions[0]* dimensions[1]* dimensions[2];
   av.initVolumeAveragingValues();

   for (int t = startTimeStep; t <= numberOfTimeSteps; t+=timeStep)
   {
      av.createMQMatrix(pathIn + "/mq/mq"+UbSystem::toString(t)+".pvtu", deltax, geo_origin);
      av.volumeAveragingWithMPI(real_l, deltax);
      av.sumOfVolumeAveragingValues();
      av.writeVolumeAveragingValuesToBinaryFiles(pathOut + "/va/vav/vav", t);
   }

   av.initMeanVolumeAveragingValues();
   av.meanOfVolumeAveragingValues(numberOfTimeSteps);
   av.writeMeanVolumeAveragingValuesToBinaryFiles(pathOut + "/va/mean/mean");
   av.initFluctuationsofVolumeAveragingValues();
   av.initSumOfFluctuations();
   av.initMeanOfFluctuations();
   av.initStresses();
   av.initSumOfStresses();
   av.initMeanOfStresses();
   av.initPlanarAveragingMQ();

   for (int t = startTimeStep; t <= numberOfTimeSteps; t+=timeStep)
   {
      av.readVolumeAveragingValuesFromBinaryFiles(pathOut + "/va/vav/vav", t);
      av.fluctuationsOfVolumeAveragingValue();
      av.writeFluctuationsToBinaryFiles(pathOut + "/va/av/Fluc", t);
      av.writeStressesToBinaryFiles(pathOut + "/va/av/Stresses", t);
      av.sumOfFluctuations();
      av.SumOfStresses();
   }

   av.meanOfFluctuations(numberOfTimeSteps);
   av.MeanOfStresses(numberOfTimeSteps);
   av.PlanarAveragingMQ(dimensions);
   av.WriteToCSV(pathOut + "/va/av", geo_origin[2], deltax);
}