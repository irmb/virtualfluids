#ifndef MEASURE_POINT_WRITER_H
#define MEASURE_POINT_WRITER_H

//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
// #include <math.h>

//#include <cmath>

//#include "Calculation/Calculation.h"
//#include "lbm/constants/D3Q27.h"
#include <numeric>
#include "basics/utilities/UbFileOutputASCII.h"
#include "Parameter/Parameter.h"

class MeasurePointWriter
{
public:
    MeasurePointWriter(){}
    ~MeasurePointWriter(){}

    static void writeMeasurePoints(Parameter* para, int level, int index, int t)
    {
        std::ostringstream convert;   // stream used for the conversion
        convert << t;      // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        UbFileOutputASCII out(para->getFName()+"_MeasurePoint_"+para->getParH(level)->MP[index].name+"_"+st+".dat");

        out.writeString("Level:");
        out.writeInteger(level);
        out.writeLine();
        out.writeString("Vx  Vy  Vz  Rho");
        out.writeLine();
        int numberNodes = (int)para->getParH(level)->MP[index].Rho.size();
        out.writeInteger(numberNodes);
        out.writeLine();
        for(int u=0; u<numberNodes; u++)
        {
            out.writeFloat((float)(para->getParH(level)->MP[index].Vx[u] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->MP[index].Vy[u] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->MP[index].Vz[u] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->MP[index].Rho[u] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeLine();
        }
        out.writeLine();
    }

    static void writeTestAcousticXY(Parameter* para, int level, int t)
    {
        std::ostringstream convert;   // stream used for the conversion
        convert << t;      // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        std::ostringstream convertLevel;   // stream used for the conversion
        convertLevel << level;      // insert the textual representation of 'Number' in the characters in the stream
        std::string sLevel = convertLevel.str();

        UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_XY.dat");

        int numberNodes = (int)para->getParH(level)->numberOfNodes;
        
        real deltaX = 1.0f / pow(2, level);
        real halfDx = deltaX / 2.0f;
        real middleOfTheGrid = (para->getMaxCoordZ()[0] + para->getMinCoordZ()[0]) / 2.f;
        //cout << "deltax: " << deltaX << ", halfDx: " << halfDx << ", middle of the grid: " << middleOfTheGrid << endl;

        for (int u = 0; u < numberNodes; u++)
        {
            if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
                ((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateZ[u]) &&
                ((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateZ[u]) )
            {
                out.writeFloat((float)(para->getParH(level)->rho[u]));
                out.writeFloat((float)(para->getParH(level)->pressure[u]));
                out.writeFloat((float)(para->getParH(level)->coordinateX[u]));
                out.writeFloat((float)(para->getParH(level)->coordinateY[u]));
                out.writeLine();
            }
        }
    }

    static void writeTestAcousticYZ(Parameter* para, int level, int t)
    {
        std::ostringstream convert;   // stream used for the conversion
        convert << t;      // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        std::ostringstream convertLevel;   // stream used for the conversion
        convertLevel << level;      // insert the textual representation of 'Number' in the characters in the stream
        std::string sLevel = convertLevel.str();

        UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_YZ.dat");

        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = 1.0f / pow(2, level);
        real halfDx = deltaX / 2.0f;
        real middleOfTheGrid = (para->getMaxCoordX()[0] + para->getMinCoordX()[0]) / 2.f;

        for (int u = 0; u < numberNodes; u++)
        {
            if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
                ((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateX[u]) &&
                ((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateX[u]))
            {
                out.writeFloat((float)(para->getParH(level)->rho[u]));
                out.writeFloat((float)(para->getParH(level)->pressure[u]));
                out.writeFloat((float)(para->getParH(level)->coordinateY[u]));
                out.writeFloat((float)(para->getParH(level)->coordinateZ[u]));
                out.writeLine();
            }
        }
    }

    static void writeTestAcousticXZ(Parameter* para, int level, int t)
    {
        std::ostringstream convert;   // stream used for the conversion
        convert << t;      // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        std::ostringstream convertLevel;   // stream used for the conversion
        convertLevel << level;      // insert the textual representation of 'Number' in the characters in the stream
        std::string sLevel = convertLevel.str();

        UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_XZ.dat");

        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = 1.0f / pow(2, level);
        real halfDx = deltaX / 2.0f;
        real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

        for (int u = 0; u < numberNodes; u++)
        {
            if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
                ((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
                ((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u]))
            {
                out.writeFloat((float)(para->getParH(level)->rho[u]));
                out.writeFloat((float)(para->getParH(level)->pressure[u]));
                out.writeFloat((float)(para->getParH(level)->coordinateX[u]));
                out.writeFloat((float)(para->getParH(level)->coordinateZ[u]));
                out.writeLine();
            }
        }
    }

    static void calcAndWriteMeanAndFluctuations(Parameter* para, int level, int t, unsigned int startOfCalculation)
    {
        //calc
        int numberNodes = (int)para->getParH(level)->MP.size();
        std::vector<double> uMean(numberNodes);
        std::vector<double> vMean(numberNodes);
        std::vector<double> wMean(numberNodes);
        std::vector<double> uuMean(numberNodes);
        std::vector<double> vvMean(numberNodes);
        std::vector<double> wwMean(numberNodes);
        std::vector<double> uvMean(numberNodes);
        std::vector<double> uwMean(numberNodes);
        std::vector<double> vwMean(numberNodes);

        int numberSteps = (int)para->getParH(level)->MP[0].Vx.size();

        std::vector<double> uuFluct(t);
        std::vector<double> vvFluct(t);
        std::vector<double> wwFluct(t);
        std::vector<double> uvFluct(t);
        std::vector<double> uwFluct(t);
        std::vector<double> vwFluct(t);

        for (int index = 0; index < numberNodes; index++)
        {
            //uMean[index] = std::accumulate(para->getParH(level)->MP[index].Vx.begin() + (startOfCalculation - 1), para->getParH(level)->MP[index].Vx.end(), 0.0) / double(t - startOfCalculation);
            //vMean[index] = std::accumulate(para->getParH(level)->MP[index].Vy.begin() + (startOfCalculation - 1), para->getParH(level)->MP[index].Vy.end(), 0.0) / double(t - startOfCalculation);
            //wMean[index] = std::accumulate(para->getParH(level)->MP[index].Vz.begin() + (startOfCalculation - 1), para->getParH(level)->MP[index].Vz.end(), 0.0) / double(t - startOfCalculation);

            uMean[index] = 0.0;
            vMean[index] = 0.0;
            wMean[index] = 0.0;

            for (int step = startOfCalculation; step < t; step++)
            {
                uMean[index] += para->getParH(level)->MP[index].Vx[step];
                vMean[index] += para->getParH(level)->MP[index].Vy[step];
                wMean[index] += para->getParH(level)->MP[index].Vz[step];
            }

            uMean[index] /= double(t - startOfCalculation);
            vMean[index] /= double(t - startOfCalculation);
            wMean[index] /= double(t - startOfCalculation);

            for (int step = 0; step < t; step++)
            {
                uuFluct[step] = (para->getParH(level)->MP[index].Vx[step] - uMean[index]) * (para->getParH(level)->MP[index].Vx[step] - uMean[index]);
                vvFluct[step] = (para->getParH(level)->MP[index].Vy[step] - vMean[index]) * (para->getParH(level)->MP[index].Vy[step] - vMean[index]);
                wwFluct[step] = (para->getParH(level)->MP[index].Vz[step] - wMean[index]) * (para->getParH(level)->MP[index].Vz[step] - wMean[index]);

                uvFluct[step] = (para->getParH(level)->MP[index].Vx[step] - uMean[index]) * (para->getParH(level)->MP[index].Vy[step] - vMean[index]);
                uwFluct[step] = (para->getParH(level)->MP[index].Vx[step] - uMean[index]) * (para->getParH(level)->MP[index].Vz[step] - wMean[index]);
                vwFluct[step] = (para->getParH(level)->MP[index].Vy[step] - vMean[index]) * (para->getParH(level)->MP[index].Vz[step] - wMean[index]);
            }

            //uuMean[index] = std::accumulate(uuFluct.begin() + (startOfCalculation - 1), uuFluct.end(), 0.0) / double(t - startOfCalculation);
            //vvMean[index] = std::accumulate(vvFluct.begin() + (startOfCalculation - 1), vvFluct.end(), 0.0) / double(t - startOfCalculation);
            //wwMean[index] = std::accumulate(wwFluct.begin() + (startOfCalculation - 1), wwFluct.end(), 0.0) / double(t - startOfCalculation);

            //uvMean[index] = std::accumulate(uvFluct.begin() + (startOfCalculation - 1), uvFluct.end(), 0.0) / double(t - startOfCalculation);
            //uwMean[index] = std::accumulate(uwFluct.begin() + (startOfCalculation - 1), uwFluct.end(), 0.0) / double(t - startOfCalculation);
            //vwMean[index] = std::accumulate(vwFluct.begin() + (startOfCalculation - 1), vwFluct.end(), 0.0) / double(t - startOfCalculation);

            uuMean[index] = 0.0;
            vvMean[index] = 0.0;
            wwMean[index] = 0.0;

            uvMean[index] = 0.0;
            uwMean[index] = 0.0;
            vwMean[index] = 0.0;

            for (int step = startOfCalculation; step < t; step++)
            {
                uuMean[index] += uuFluct[step];
                vvMean[index] += vvFluct[step];
                wwMean[index] += wwFluct[step];

                uvMean[index] += uvFluct[step];
                uwMean[index] += uwFluct[step];
                vwMean[index] += vwFluct[step];
            }

            uuMean[index] /= double(t - startOfCalculation);
            vvMean[index] /= double(t - startOfCalculation);
            wwMean[index] /= double(t - startOfCalculation);

            uvMean[index] /= double(t - startOfCalculation);
            uwMean[index] /= double(t - startOfCalculation);
            vwMean[index] /= double(t - startOfCalculation);

        }

        std::cout << "level: " << level << ", t: " << t << ", startOfCalculation= " << startOfCalculation << ", numberSteps: " << numberSteps << std::endl;

        //write
        UbFileOutputASCII out(para->getFName() + "_MeanAndFluct_" + std::to_string(level) + "_" + std::to_string(t) + ".dat");

        out.writeString("Level:");
        out.writeInteger(level);
        out.writeLine();
        out.writeString("PointName uMean  vMean  wMean  uuMean  vvMean  wwMean  uvMean  uwMean  vwMean");
        out.writeLine();
        out.writeInteger(numberNodes);
        out.writeLine();
        for (int index = 0; index < numberNodes; index++)
        {
            out.writeString(para->getParH(level)->MP[index].name);
            out.writeFloat((float)(uMean[index]  * para->getVelocityRatio()));
            out.writeFloat((float)(vMean[index]  * para->getVelocityRatio()));
            out.writeFloat((float)(wMean[index]  * para->getVelocityRatio()));
            out.writeFloat((float)(uuMean[index] * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeFloat((float)(vvMean[index] * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeFloat((float)(wwMean[index] * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeFloat((float)(uvMean[index] * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeFloat((float)(uwMean[index] * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeFloat((float)(vwMean[index] * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeLine();
        }
        out.writeLine();
    }

    static void writeNodes(SPtr<Parameter> para, int level)
    {
        std::ostringstream convert;   // stream used for the conversion
        std::string st = convert.str();
        UbFileOutputASCII out(para->getFName() + "_Nodes_" + std::to_string(level) + ".dat");

        out.writeString("Level:");
        out.writeInteger(level);
        out.writeLine();
        out.writeString("Index  Nx  Ny  Ny");
        out.writeLine();
        int numberNodes = (int)para->getParH(level)->numberOfNodes;
        out.writeInteger(numberNodes);
        out.writeLine();
        for (int u = 0; u < numberNodes; u++)
        {
            out.writeInteger((int)(u));
            out.writeInteger((int)(para->getParH(level)->neighborX[u]));
            out.writeInteger((int)(para->getParH(level)->neighborY[u]));
            out.writeInteger((int)(para->getParH(level)->neighborZ[u]));
            out.writeLine();
        }
        out.writeLine();
    }
protected:

private:
};
#endif
