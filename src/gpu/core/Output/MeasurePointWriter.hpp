//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================
#ifndef MEASURE_POINT_WRITER_H
#define MEASURE_POINT_WRITER_H

#include <numeric>

#include <basics/utilities/UbFileOutputASCII.h>

#include "Parameter/Parameter.h"

class MeasurePointWriter
{
public:
    static void writeMeasurePoints(Parameter* para, int level, int index, int t)
    {
        std::ostringstream convert; // stream used for the conversion
        convert << t;               // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        UbFileOutputASCII out(para->getFName() + "_MeasurePoint_" + para->getParH(level)->MeasurePointVector[index].name + "_" + st +
                              ".dat");

        out.writeString("Level:");
        out.writeInteger(level);
        out.writeLine();
        out.writeString("Vx  Vy  Vz  Rho");
        out.writeLine();
        int numberNodes = (int)para->getParH(level)->MeasurePointVector[index].Rho.size();
        out.writeInteger(numberNodes);
        out.writeLine();
        for (int u = 0; u < numberNodes; u++) {
            out.writeFloat((float)(para->getParH(level)->MeasurePointVector[index].Vx[u] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->MeasurePointVector[index].Vy[u] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->MeasurePointVector[index].Vz[u] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->MeasurePointVector[index].Rho[u] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio()));
            out.writeLine();
        }
        out.writeLine();
    }

    static void writeTestAcousticXY(Parameter* para, int level, int t)
    {
        std::ostringstream convert; // stream used for the conversion
        convert << t;               // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        std::ostringstream convertLevel; // stream used for the conversion
        convertLevel << level;           // insert the textual representation of 'Number' in the characters in the stream
        std::string sLevel = convertLevel.str();

        UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_XY.dat");

        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = vf::basics::constant::c1o1 / pow(2, level);
        real halfDx = deltaX / 2.0f;
        real middleOfTheGrid = (para->getMaxCoordZ()[0] + para->getMinCoordZ()[0]) / vf::basics::constant::c2o1;

        for (int u = 0; u < numberNodes; u++) {
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
        std::ostringstream convert; // stream used for the conversion
        convert << t;               // insert the textual representation of 'Number' in the characters in the stream
        std::string st = convert.str();
        std::ostringstream convertLevel; // stream used for the conversion
        convertLevel << level;           // insert the textual representation of 'Number' in the characters in the stream
        std::string sLevel = convertLevel.str();

        UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_YZ.dat");

        int numberNodes = (int)para->getParH(level)->numberOfNodes;

        real deltaX = vf::basics::constant::c1o1 / pow(2, level);
        real halfDx = deltaX / vf::basics::constant::c2o1;
        real middleOfTheGrid = (para->getMaxCoordX()[0] + para->getMinCoordX()[0]) / vf::basics::constant::c2o1;

        for (int u = 0; u < numberNodes; u++) {
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
        int numberNodes = (int)para->getParH(level)->MeasurePointVector.size();
        std::vector<double> uMean(numberNodes);
        std::vector<double> vMean(numberNodes);
        std::vector<double> wMean(numberNodes);
        std::vector<double> uuMean(numberNodes);
        std::vector<double> vvMean(numberNodes);
        std::vector<double> wwMean(numberNodes);
        std::vector<double> uvMean(numberNodes);
        std::vector<double> uwMean(numberNodes);
        std::vector<double> vwMean(numberNodes);

        int numberSteps = (int)para->getParH(level)->MeasurePointVector[0].Vx.size();

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
                uMean[index] += para->getParH(level)->MeasurePointVector[index].Vx[step];
                vMean[index] += para->getParH(level)->MeasurePointVector[index].Vy[step];
                wMean[index] += para->getParH(level)->MeasurePointVector[index].Vz[step];
            }

            uMean[index] /= double(t - startOfCalculation);
            vMean[index] /= double(t - startOfCalculation);
            wMean[index] /= double(t - startOfCalculation);

            for (int step = 0; step < t; step++)
            {
                uuFluct[step] = (para->getParH(level)->MeasurePointVector[index].Vx[step] - uMean[index]) * (para->getParH(level)->MeasurePointVector[index].Vx[step] - uMean[index]);
                vvFluct[step] = (para->getParH(level)->MeasurePointVector[index].Vy[step] - vMean[index]) * (para->getParH(level)->MeasurePointVector[index].Vy[step] - vMean[index]);
                wwFluct[step] = (para->getParH(level)->MeasurePointVector[index].Vz[step] - wMean[index]) * (para->getParH(level)->MeasurePointVector[index].Vz[step] - wMean[index]);

                uvFluct[step] = (para->getParH(level)->MeasurePointVector[index].Vx[step] - uMean[index]) * (para->getParH(level)->MeasurePointVector[index].Vy[step] - vMean[index]);
                uwFluct[step] = (para->getParH(level)->MeasurePointVector[index].Vx[step] - uMean[index]) * (para->getParH(level)->MeasurePointVector[index].Vz[step] - wMean[index]);
                vwFluct[step] = (para->getParH(level)->MeasurePointVector[index].Vy[step] - vMean[index]) * (para->getParH(level)->MeasurePointVector[index].Vz[step] - wMean[index]);
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
            out.writeString(para->getParH(level)->MeasurePointVector[index].name);
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
};

#endif
