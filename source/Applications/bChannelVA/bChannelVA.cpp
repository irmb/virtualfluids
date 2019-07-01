#include <iostream>
#include <string>
#include "VirtualFluids.h"

#include "Averaging.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   try
   {
      SPtr<Communicator> comm = MPICommunicator::getInstance();
      int myid = comm->getProcessID();

      //Pheonix
      //double deltaX = 1;
      //double halfDeltaX = deltaX / 2.0;
      //std::array<int, 3> dimensions = { 600 / (int)deltaX, 400 / (int)deltaX, 400 / (int)deltaX };
      //std::array<double, 3> geo_origin = { halfDeltaX, halfDeltaX, halfDeltaX };
      //std::array<double, 3> geo_spacing = { 1,1,1 };
      //std::array<int, 6> geo_extent = { 0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1 };
      //double real_l = 40;
      //double l = 40;

      //int startTimeStep = 600000;
      //int timeStep = 10000;
      //int numberOfTimeSteps = 610000; //1200000;
      //int numberOfSamples = numberOfTimeSteps / startTimeStep;
      //int numberOfGridPoints = dimensions[0] * dimensions[1] * dimensions[2];


      //Bombadil
      string pathIn = "d:/temp/BreugemChannelAnisotrop2";
      string pathOut = "d:/temp/BreugemChannelAnisotrop2";

      double deltaX = 10;
      double halfDeltaX = deltaX / 2.0;
      std::array<int, 3> dimensions = { 600 / (int)deltaX, 400 / (int)deltaX, 400 / (int)deltaX };
      std::array<double, 3> geo_origin = { halfDeltaX, halfDeltaX, halfDeltaX };
      std::array<double, 3> geo_spacing = { 10,10,10 };
      std::array<int, 6> geo_extent = { 0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1 };
      double real_l = 20;
      double l = 20;

      int startTimeStep = 60000;
      int timeStep = 1000;
      int numberOfTimeSteps = 65000;
      int numberOfSamples = (numberOfTimeSteps - startTimeStep) / timeStep + 1;
      int numberOfGridPoints = dimensions[0] * dimensions[1] * dimensions[2];

      Averaging av;

      av.setDimensions(dimensions);
      av.setExtent(geo_extent);
      av.setOrigin(geo_origin);
      av.setSpacing(geo_spacing);
      av.setDeltaX(deltaX);

      //read geo matrix
      av.createGeoMatrix(pathIn + "/bc/bc0.pvtu");
      if (myid == 0) av.writeGeoMatrixToBinaryFiles(pathOut + "/va/geo/geomatrix.bin");
      av.readGeoMatrixFromBinaryFiles(pathOut + "/va/geo/geomatrix.bin");

      //read mq matrix
      for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
      {
         av.createMQMatrix(pathIn + "/mq/mq" + UbSystem::toString(t) + ".pvtu");
         av.writeMqMatrixToBinaryFiles(pathOut + "/va/mq/mq", t);
      }

      //compute mq values
      av.initMeanMqValues();
      for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
      {
         av.readMqMatrixFromBinaryFiles(pathOut + "/va/mq/mq", t);
         av.sumMqValues();
      }
      av.computeMeanMqValues(numberOfSamples);
      av.writeMeanMqValuesToBinaryFiles(pathOut + "/va/mean/mean");

      //compute volume averaging of Reynolds averaged MQ values
      av.volumeAveragingOfMeanMqValuesWithMPI(l);
      av.writeVaMeanMqValuesToBinaryFiles(pathOut + "/va/vaMean/vaMean");

      //compute fluctuations
      av.readMeanMqValuesFromBinaryFiles(pathOut + "/va/mean/mean");
      av.initFluctuationsOfMqValues();
      for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
      {
         av.readMqMatrixFromBinaryFiles(pathOut + "/va/mq/mq", t);
         av.computeFluctuationsOfMqValues();
         av.writeFluctuationsOfMqValuesToBinaryFiles(pathOut + "/va/fluc/fluc", t);
      }

      //compute volume averaged fluctuations
      av.initMeanOfVolumeAveragedValues();
      for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
      {
         av.readFluctuationsOfMqValuesFromBinaryFiles(pathOut + "/va/fluc/fluc", t);
         av.volumeAveragingOfFluctuationsWithMPI(l);
         av.writeVaFluctuationsToBinaryFiles(pathOut + "/va/vaFluc/vaFluc", t);
         av.sumVolumeAveragedValues();
      }
      av.computeVolumeAveragedValues(numberOfSamples);
      av.writeVolumeAveragedValuesToBinaryFiles(pathOut + "/va/values/val");

      //planar averaging
      av.initPlanarAveraging();
      av.planarAveraging();
      av.writeToCSV(pathOut + "/va/planar/planar", geo_origin[2], deltaX);
   }
   catch (const std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (std::string& s)
   {
      cerr << s << endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
   }

}