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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "AnalyticalResults2DToVTKWriterImp.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"

#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"

#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include <mpi.h>
#include "core/Output/FilePartCalculator.h"

std::shared_ptr<AnalyticalResults2DToVTKWriterImp>
AnalyticalResults2DToVTKWriterImp::getInstance(bool writeAnalyticalResults)
{
    static std::shared_ptr<AnalyticalResults2DToVTKWriterImp> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<AnalyticalResults2DToVTKWriterImp>(
            new AnalyticalResults2DToVTKWriterImp(writeAnalyticalResults));
    return uniqueInstance;
}

AnalyticalResults2DToVTKWriterImp::AnalyticalResults2DToVTKWriterImp(bool writeAnalyticalResults)
    : writeAnalyticalResults(writeAnalyticalResults)
{
}

void AnalyticalResults2DToVTKWriterImp::writeAnalyticalResult(std::shared_ptr<Parameter> para,
                                                              std::shared_ptr<AnalyticalResults> analyticalResult)
{
    if (writeAnalyticalResults) {
        std::cout << "Write Analytical Result To VTK-Files" << std::endl;
        for (int level = para->getCoarse(); level <= para->getFine(); level++) {
#pragma omp parallel for
            for (int timeStep = 0; timeStep < analyticalResult->getNumberOfTimeSteps(); timeStep++) {
                const unsigned int numberOfParts = FilePartCalculator::calculateNumberOfParts(para->getParH(level)->numberOfNodes);
                std::vector<std::string> fname;
                unsigned int time =
                    analyticalResult->getTimeSteps().at(timeStep) * analyticalResult->getTimeStepLength();
                for (int j = 1; j <= numberOfParts; j++) {
                    std::string filePath = para->getFName();
                    filePath.resize(filePath.size() - 5);
                    fname.push_back(filePath + "AnalyticalResult/Analytical_cells_bin_lev_" +
                                    StringUtil::toString<int>(level) + "_ID_" +
                                    StringUtil::toString<int>(para->getMyProcessID()) + "_Part_" +
                                    StringUtil::toString<int>(j) + "_t_" + StringUtil::toString<int>(time) + ".vtk");
                }
                std::cout << "\t Write TimeStep=" << timeStep << " t=" << time << "...";
                writeTimeStep(para, analyticalResult, level, fname, timeStep);
                std::cout << "done." << std::endl;
            }
        }
        std::cout << std::endl;
    }
}

void AnalyticalResults2DToVTKWriterImp::writeTimeStep(std::shared_ptr<Parameter> para,
                                                      std::shared_ptr<AnalyticalResults> analyticalResult, int level,
                                                      std::vector<std::string> &fname, int timeStep)
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    std::vector<std::vector<double>> nodedata(nodedatanames.size());

    maxX = para->getGridX().at(level);
    maxY = para->getGridY().at(level);
    maxZ = para->getGridZ().at(level);

    std::vector<double> press = analyticalResult->getPress()[timeStep];
    std::vector<double> rho = analyticalResult->getRho()[timeStep];
    std::vector<double> vx = analyticalResult->getVx()[timeStep];
    std::vector<double> vy = analyticalResult->getVy()[timeStep];
    std::vector<double> vz = analyticalResult->getVz()[timeStep];

    for (unsigned int part = 0; part < fname.size(); part++) {
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes, part);

        //////////////////////////////////////////////////////////////////////////
        startpos = FilePartCalculator::calculateStartingPostionOfPart(part);
        endpos = startpos + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        nodedata[0].resize(sizeOfNodes);
        nodedata[1].resize(sizeOfNodes);
        nodedata[2].resize(sizeOfNodes);
        nodedata[3].resize(sizeOfNodes);
        nodedata[4].resize(sizeOfNodes);
        nodedata[5].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        for (unsigned int pos = startpos; pos < endpos; pos++) {
            std::cout << "BEGIN POS: " << pos << std::endl;
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID) {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsAreFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));

                int numberInResults = CoordResults2DTo1D(x1 - 1.0, x3 - 1.0);
                nodedata[0][dn1] = press[numberInResults];
                nodedata[1][dn1] = rho[numberInResults];
                nodedata[2][dn1] = vx[numberInResults];
                nodedata[3][dn1] = vy[numberInResults];
                nodedata[4][dn1] = vz[numberInResults];
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                std::cout << "NeighborIndex1 " << number1 << std::endl <<
                "NeighborIndex2 " << number2 << std::endl <<
                "NeighborIndex3 " << number3 << std::endl <<
                "NeighborIndex4 " << number4 << std::endl <<
                "NeighborIndex5 " << number5 << std::endl <<
                "NeighborIndex6 " << number6 << std::endl <<
                "NeighborIndex7 " << number7 << std::endl <<
                "NeighborIndex8 " << number8 << std::endl;
                //////////////////////////////////////////////////////////////////////////
                auto neighbor1 = para->getParH(level)->typeOfGridNode[number2];
                auto neighbor2 = para->getParH(level)->typeOfGridNode[number3];
                auto neighbor3 = para->getParH(level)->typeOfGridNode[number4];
                auto neighbor4 = para->getParH(level)->typeOfGridNode[number5]; //breaks!
                auto neighbor5 = para->getParH(level)->typeOfGridNode[number6];
                auto neighbor6 = para->getParH(level)->typeOfGridNode[number7];
                auto neighbor7 = para->getParH(level)->typeOfGridNode[number8];
                std::cout << "Neighbor1 " << neighbor1 << std::endl <<
                "Neighbor2 " << neighbor2 << std::endl <<
                "Neighbor3 " << neighbor3 << std::endl <<
                "Neighbor4 " << neighbor4 << std::endl <<
                "Neighbor5 " << neighbor5 << std::endl <<
                "Neighbor6 " << neighbor6 << std::endl <<
                "Neighbor7 " << neighbor7 << std::endl;

                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID)
                    neighborsAreFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos || number3 > endpos || number4 > endpos || number5 > endpos || number6 > endpos ||
                    number7 > endpos || number8 > endpos)
                    neighborsAreFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsAreFluid)
                    cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                std::cout << "END POS: " << pos << std::endl;
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
    }
}

int AnalyticalResults2DToVTKWriterImp::CoordResults2DTo1D(int x, int z)
{
    return z * (maxX - 1) + x;
}

//! \}
