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
        string pathIn = "d:/temp/BreugemChannelAnisotrop";
        string pathOut = "d:/temp/BreugemChannelAnisotrop";

        double deltaX = 10;
        double halfDeltaX = deltaX / 2.0;
        std::array<int, 3> dimensions = { 600/(int)deltaX, 400/(int)deltaX, 400/(int)deltaX };
        std::array<double, 3> geo_origin = { halfDeltaX, halfDeltaX, halfDeltaX };
        std::array<double, 3> geo_spacing = { 10,10,10 };
        std::array<int, 6> geo_extent = { 0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1 };
        double real_l = 40;
        double l = 40;

        int startTimeStep = 100;
        int timeStep = 100;
        int numberOfTimeSteps = 500;
        int numberOfSamples = (numberOfTimeSteps - startTimeStep)/timeStep + 1;
        int numberOfGridPoints = dimensions[0] * dimensions[1] * dimensions[2];

        Averaging av;
        av.setDimensions(dimensions);
        av.setExtent(geo_extent);
        av.setOrigin(geo_origin);
        av.setSpacing(geo_spacing);
        av.setDeltaX(deltaX);

        av.createGeoMatrix(pathIn + "/bc/bc0.pvtu");
        if (myid == 0) av.writeGeoMatrixToBinaryFiles(pathOut + "/va/geomatrix.bin");
        //av.readGeoMatrixFromBinaryFiles(pathOut + "/va/geomatrix.bin");
        av.initVolumeAveragingValues();

        for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
        {
            av.createMQMatrix(pathIn + "/mq/mq" + UbSystem::toString(t) + ".pvtu");
            if (myid == 0) av.writeMqMatrixToBinaryFiles(pathOut + "/va/mq/mq" + UbSystem::toString(t) + "/mq", t);
            av.volumeAveragingWithMPI(real_l);
            if (myid == 0)
            {
               av.sumOfVolumeAveragingValues();
               av.writeVolumeAveragingValuesToBinaryFiles(pathOut + "/va/vav/vav", t);
            }
        }

        //av.writeMqMatrixToImageFile(pathOut + "/va/vav/mq");
        //if (myid == 0) av.readVolumeAveragingValuesFromBinaryFiles(pathOut + "/va/vav/vav", 60000);
        //if (myid == 0) av.writeVaMatrixToImageFile(pathOut + "/va/vav/vav");

        if (myid == 0)
        {
           av.initMeanVolumeAveragingValues();
           av.meanOfVolumeAveragingValues(numberOfSamples);
           av.writeMeanVolumeAveragingValuesToBinaryFiles(pathOut + "/va/mean/mean");
           av.writeMeanMatrixToImageFile(pathOut + "/va/vav/vavMean");
        }

        //for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
        //{
        //   av.readVolumeAveragingValuesFromBinaryFiles(pathOut + "/va/vav/vav", t);
        //   av.sumOfVolumeAveragingValues();
        //}

        //av.initMeanVolumeAveragingValues();
        //av.meanOfVolumeAveragingValues(numberOfSamples);
        //av.writeMeanVolumeAveragingValuesToBinaryFiles(pathOut + "/va/mean/mean");
        //av.writeMeanMatrixToImageFile(pathOut + "/va/vav/vavMean");

        av.initFluctuations();
        av.initSumOfVaFluctuations();
        av.initMeanOfVaFluctuations();

        av.initStresses();
        av.initSumOfVaStresses();
        av.initMeanOfVaStresses();
        av.initPlanarAveraging();

        av.initVolumeAveragingFluctStressValues();

        for (int t = startTimeStep; t <= numberOfTimeSteps; t += timeStep)
        {
            //av.readVolumeAveragingValuesFromBinaryFiles(pathOut + "/va/vav/vav", t);
            //av.readMqMatrixFromBinaryFiles(pathOut + "/va/mq/mq" + UbSystem::toString(t) + "/mq", t);
            av.fluctuationsStress();
            av.writeFluctuationsToImageFile(pathOut + "/va/vav/vaFluct");

            av.volumeAveragingFluctStressWithMPI(real_l);
            av.writeVaFluctuationsToBinaryFiles(pathOut + "/va/av/Fluc", t);
            av.writeVaStressesToBinaryFiles(pathOut + "/va/av/Stresses", t);
            av.writeVaFluctuationsToImageFile(pathOut + "/va/vav/vaVAFluct");

            av.sumOfVaFluctuations();
            av.sumOfVaStresses();
        }

        av.meanOfVaFluctuations(numberOfSamples);
        av.meanOfVaStresses(numberOfSamples);
        av.writeMeanOfVaFluctuationsToImageFile(pathOut + "/va/vav/vaMeanFluct");

        av.planarAveraging();
        av.writeToCSV(pathOut + "/va/av", geo_origin[2], deltaX);
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