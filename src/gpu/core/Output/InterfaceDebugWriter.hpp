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
#ifndef INTERFACEDEBUG_HPP
#define INTERFACEDEBUG_HPP

#include <fstream>
#include <sstream>
#include <cmath>

#include "StringUtilities/StringUtil.h"
#include "lbm/constants/D3Q27.h"
#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"
#include <basics/writer/WbWriterVtkXmlBinary.h>

namespace InterfaceDebugWriter
{

void writeGridInterfaceLines(Parameter *para, int level, const uint *coarse, const uint *fine, uint numberOfNodes,
                             const std::string &name)
{
    std::vector<UbTupleFloat3> nodes(numberOfNodes * 2);
    std::vector<UbTupleInt2> cells(numberOfNodes);

    int actualNodeNumber = 0;
    for (uint u = 0; u < numberOfNodes; u++) {
        const int posCoarse   = coarse[u];
        const double x1Coarse = para->getParH(level)->coordinateX[posCoarse];
        const double x2Coarse = para->getParH(level)->coordinateY[posCoarse];
        const double x3Coarse = para->getParH(level)->coordinateZ[posCoarse];

        const int posFine   = fine[u];
        const double x1Fine = para->getParH(level + 1)->coordinateX[posFine];
        const double x2Fine = para->getParH(level + 1)->coordinateY[posFine];
        const double x3Fine = para->getParH(level + 1)->coordinateZ[posFine];

        nodes[actualNodeNumber++] = makeUbTuple(float(x1Coarse), float(x2Coarse), float(x3Coarse));
        nodes[actualNodeNumber++] = makeUbTuple(float(x1Fine), float(x2Fine), float(x3Fine));

        cells[u] = makeUbTuple(actualNodeNumber - 2, actualNodeNumber - 1);
    }
    WbWriterVtkXmlBinary::getInstance()->writeLines(name, nodes, cells);
}

void writeInterfaceLinesDebugCF(Parameter *para)
{
    for (int level = 0; level < para->getMaxLevel(); level++) {
        const std::string fileName = para->getFName() + "_" + StringUtil::toString<int>(level) + "_OffDebugCF.vtk";
        writeGridInterfaceLines(para, level, para->getParH(level)->coarseToFine.coarseCellIndices, para->getParH(level)->coarseToFine.fineCellIndices, para->getParH(level)->coarseToFine.numberOfCells, fileName);
    }
}

void writeInterfaceLinesDebugFC(Parameter *para)
{
    for (int level = 0; level < para->getMaxLevel(); level++) {
        const std::string fileName = para->getFName() + "_" + StringUtil::toString<int>(level) + "_OffDebugFC.vtk";
        writeGridInterfaceLines(para, level, para->getParH(level)->fineToCoarse.coarseCellIndices, para->getParH(level)->fineToCoarse.fineCellIndices, para->getParH(level)->fineToCoarse.numberOfCells, fileName);
    }
}

//////////////////////////////////////////////////////////////////////////
void writeGridInterfaceLinesNeighbors(Parameter *para, int level, const uint *interfaceIndices, uint numberOfNodes,
                                      const std::string &name)
{
    std::vector<UbTupleFloat3> nodes(numberOfNodes * 2);
    std::vector<UbTupleInt2> cells(numberOfNodes);

    int actualNodeNumber = 0;
    for (uint u = 0; u < numberOfNodes; u++) {
        const int pos   = interfaceIndices[u];
        const double x1 = para->getParH(level)->coordinateX[pos];
        const double x2 = para->getParH(level)->coordinateY[pos];
        const double x3 = para->getParH(level)->coordinateZ[pos];

        const double x1Neighbor = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
        const double x2Neighbor = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
        const double x3Neighbor = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];

        nodes[actualNodeNumber++] = (makeUbTuple(float(x1), float(x2), float(x3)));
        nodes[actualNodeNumber++] = (makeUbTuple(float(x1Neighbor), float(x2Neighbor), float(x3Neighbor)));

        cells[u] = makeUbTuple(actualNodeNumber - 2, actualNodeNumber - 1);
    }
    WbWriterVtkXmlBinary::getInstance()->writeLines(name, nodes, cells);
}

void writeInterfaceLinesDebugCFCneighbor(Parameter *para)
{
    for (int level = 0; level < para->getMaxLevel(); level++) {
        std::string filename = para->getFName() + "_" + StringUtil::toString<int>(level) + "_CFCneighbor.vtk";
        writeGridInterfaceLinesNeighbors(para, level, para->getParH(level)->coarseToFine.coarseCellIndices, para->getParH(level)->coarseToFine.numberOfCells,
                                         filename);
    }
}

//////////////////////////////////////////////////////////////////////////
void writeInterfaceLinesDebugCFFneighbor(Parameter *para)
{
    for (int level = 0; level < para->getMaxLevel(); level++) {
        std::string filename = para->getFName() + "_" + StringUtil::toString<int>(level) + "_CFFneighbor.vtk";
        writeGridInterfaceLinesNeighbors(para, level + 1, para->getParH(level)->coarseToFine.fineCellIndices, para->getParH(level)->coarseToFine.numberOfCells, filename);
    }
}

//////////////////////////////////////////////////////////////////////////
void writeInterfaceLinesDebugFCCneighbor(Parameter *para)
{
    for (int level = 0; level < para->getMaxLevel(); level++) {
        std::string filename = para->getFName() + "_" + StringUtil::toString<int>(level) + "_FCCneighbor.vtk";
        writeGridInterfaceLinesNeighbors(para, level, para->getParH(level)->fineToCoarse.coarseCellIndices, para->getParH(level)->fineToCoarse.numberOfCells,
                                         filename);
    }
}

//////////////////////////////////////////////////////////////////////////
void writeInterfaceLinesDebugFCFneighbor(Parameter *para)
{
    for (int level = 0; level < para->getMaxLevel(); level++) {
        std::string filename = para->getFName() + "_" + StringUtil::toString<int>(level) + "_FCFneighbor.vtk";
        writeGridInterfaceLinesNeighbors(para, level + 1, para->getParH(level)->fineToCoarse.fineCellIndices, para->getParH(level)->fineToCoarse.numberOfCells, filename);
    }
}

//////////////////////////////////////////////////////////////////////////
void writeInterfaceLinesDebugOff(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<UbTupleInt2> cellsVec;
    int nodeNumberVec = 0;

    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->coarseToFine.numberOfCells;
    }
    nodesVec.resize(nodeNumberVec * 8);
    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
            double xoff = para->getParH(level)->neighborCoarseToFine.x[u];
            double yoff = para->getParH(level)->neighborCoarseToFine.y[u];
            double zoff = para->getParH(level)->neighborCoarseToFine.z[u];

            int posFine = para->getParH(level)->coarseToFine.fineCellIndices[u];

            double x1Fine = para->getParH(level + 1)->coordinateX[posFine];
            double x2Fine = para->getParH(level + 1)->coordinateY[posFine];
            double x3Fine = para->getParH(level + 1)->coordinateZ[posFine];

            double x1 = x1Fine + xoff;
            double x2 = x2Fine + yoff;
            double x3 = x3Fine + zoff;

            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1Fine), (float)(x2Fine), (float)(x3Fine)));

            cellsVec.push_back(makeUbTuple(nodeCount - 2, nodeCount - 1));
        }
        std::string filenameVec = para->getFName() + "_" + StringUtil::toString<int>(level) + "_OffDebugCF_Offs.vtk";
        WbWriterVtkXmlBinary::getInstance()->writeLines(filenameVec, nodesVec, cellsVec);
        cellsVec.clear();
        nodesVec.clear();
    }
}

//////////////////////////////////////////////////////////////////////////

void writeInterfacePointsDebugCFC(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec2;
    int nodeNumberVec = 0;

    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->coarseToFine.numberOfCells;
    }
    nodesVec2.resize(nodeNumberVec * 8);
    int nodeCount2 = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
            int pos = para->getParH(level)->coarseToFine.coarseCellIndices[u];

            double x1 = para->getParH(level)->coordinateX[pos];
            double x2 = para->getParH(level)->coordinateY[pos];
            double x3 = para->getParH(level)->coordinateZ[pos];

            nodesVec2[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec2 = para->getFName() + "_" + StringUtil::toString<int>(level) + "_OffDebugPointsCF.vtk";
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec2, nodesVec2);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeBcPointsDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec2;
    int nodeNumberVec = 0;

    for (int level = 0; level <= para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->noSlipBC.numberOfBCnodes;
    }
    nodesVec2.resize(nodeNumberVec * 8);
    int nodeCount2 = 0;
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (uint u = 0; u < para->getParH(level)->noSlipBC.numberOfBCnodes; u++) {
            int pos = para->getParH(level)->noSlipBC.k[u];

            double x1 = para->getParH(level)->coordinateX[pos];
            double x2 = para->getParH(level)->coordinateY[pos];
            double x3 = para->getParH(level)->coordinateZ[pos];

            nodesVec2[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec2 = para->getFName() + "_PointsBc_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec2, nodesVec2);
    }
}

//////////////////////////////////////////////////////////////////////////

void writePressPointsDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    for (int level = 0; level <= para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->pressureBC.numberOfBCnodes;
    }
    nodesVec.resize(nodeNumberVec);
    int nodeCount2 = 0;
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (uint u = 0; u < para->getParH(level)->pressureBC.numberOfBCnodes; u++) {
            int pos = para->getParH(level)->pressureBC.k[u];

            double x1 = para->getParH(level)->coordinateX[pos];
            double x2 = para->getParH(level)->coordinateY[pos];
            double x3 = para->getParH(level)->coordinateZ[pos];

            nodesVec[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec = para->getFName() + "_PointsPress_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec, nodesVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writePressNeighborPointsDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        nodeNumberVec += (int)para->getParH(level)->pressureBC.numberOfBCnodes;
    }
    nodesVec.resize(nodeNumberVec);
    int nodeCount2 = 0;
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (uint u = 0; u < para->getParH(level)->pressureBC.numberOfBCnodes; u++) {
            int pos = para->getParH(level)->pressureBC.kN[u];

            real x1 = para->getParH(level)->coordinateX[pos];
            real x2 = para->getParH(level)->coordinateY[pos];
            real x3 = para->getParH(level)->coordinateZ[pos];

            nodesVec[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec = para->getFName() + "_PointsPressNeighbor_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec, nodesVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeNeighborXPointsDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        nodeNumberVec += (int)para->getParH(level)->numberOfNodes;
    }
    nodesVec.resize(nodeNumberVec);
    int nodeCount2 = 0;
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
            real x1 = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[index]];
            real x2 = para->getParH(level)->coordinateY[para->getParH(level)->neighborX[index]];
            real x3 = para->getParH(level)->coordinateZ[para->getParH(level)->neighborX[index]];

            nodesVec[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec = para->getFName() + "_PointsNeighborX_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec, nodesVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeNeighborXLinesDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<UbTupleInt2> cellsVec;
    int nodeNumberVec = 0;

    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->numberOfNodes;
    }
    nodesVec.resize(nodeNumberVec * 2);
    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
            real x1  = para->getParH(level)->coordinateX[index];
            real x2  = para->getParH(level)->coordinateY[index];
            real x3  = para->getParH(level)->coordinateZ[index];
            real x1N = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[index]];
            real x2N = para->getParH(level)->coordinateY[para->getParH(level)->neighborX[index]];
            real x3N = para->getParH(level)->coordinateZ[para->getParH(level)->neighborX[index]];

            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1N), (float)(x2N), (float)(x3N)));

            if (para->getParH(level)->typeOfGridNode[index] == GEO_FLUID) {
                cellsVec.push_back(makeUbTuple(nodeCount - 2, nodeCount - 1));
            }
        }
        std::string filenameVec = para->getFName() + "_" + StringUtil::toString<int>(level) + "_NeighborX_Lines.vtk";
        WbWriterVtkXmlBinary::getInstance()->writeLines(filenameVec, nodesVec, cellsVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeNeighborYPointsDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        nodeNumberVec += (int)para->getParH(level)->numberOfNodes;
    }
    nodesVec.resize(nodeNumberVec);
    int nodeCount2 = 0;
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
            real x1 = para->getParH(level)->coordinateX[para->getParH(level)->neighborY[index]];
            real x2 = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[index]];
            real x3 = para->getParH(level)->coordinateZ[para->getParH(level)->neighborY[index]];

            nodesVec[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec = para->getFName() + "_PointsNeighborY_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec, nodesVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeNeighborYLinesDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<UbTupleInt2> cellsVec;
    int nodeNumberVec = 0;

    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->numberOfNodes;
    }
    nodesVec.resize(nodeNumberVec * 2);
    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
            real x1  = para->getParH(level)->coordinateX[index];
            real x2  = para->getParH(level)->coordinateY[index];
            real x3  = para->getParH(level)->coordinateZ[index];
            real x1N = para->getParH(level)->coordinateX[para->getParH(level)->neighborY[index]];
            real x2N = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[index]];
            real x3N = para->getParH(level)->coordinateZ[para->getParH(level)->neighborY[index]];

            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1N), (float)(x2N), (float)(x3N)));

            if (para->getParH(level)->typeOfGridNode[index] == GEO_FLUID) {
                cellsVec.push_back(makeUbTuple(nodeCount - 2, nodeCount - 1));
            }
        }
        std::string filenameVec = para->getFName() + "_" + StringUtil::toString<int>(level) + "_NeighborY_Lines.vtk";
        WbWriterVtkXmlBinary::getInstance()->writeLines(filenameVec, nodesVec, cellsVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeNeighborZPointsDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        nodeNumberVec += (int)para->getParH(level)->numberOfNodes;
    }
    nodesVec.resize(nodeNumberVec);
    int nodeCount2 = 0;
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
            real x1 = para->getParH(level)->coordinateX[para->getParH(level)->neighborZ[index]];
            real x2 = para->getParH(level)->coordinateY[para->getParH(level)->neighborZ[index]];
            real x3 = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[index]];

            nodesVec[nodeCount2++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
        }
        std::string filenameVec = para->getFName() + "_PointsNeighborZ_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec, nodesVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeNeighborZLinesDebug(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<UbTupleInt2> cellsVec;
    int nodeNumberVec = 0;

    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->numberOfNodes;
    }
    nodesVec.resize(nodeNumberVec * 2);
    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++) {
            real x1  = para->getParH(level)->coordinateX[index];
            real x2  = para->getParH(level)->coordinateY[index];
            real x3  = para->getParH(level)->coordinateZ[index];
            real x1N = para->getParH(level)->coordinateX[para->getParH(level)->neighborZ[index]];
            real x2N = para->getParH(level)->coordinateY[para->getParH(level)->neighborZ[index]];
            real x3N = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[index]];

            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1N), (float)(x2N), (float)(x3N)));

            if (para->getParH(level)->typeOfGridNode[index] == GEO_FLUID) {
                cellsVec.push_back(makeUbTuple(nodeCount - 2, nodeCount - 1));
            }
        }
        std::string filenameVec = para->getFName() + "_" + StringUtil::toString<int>(level) + "_NeighborZ_Lines.vtk";
        WbWriterVtkXmlBinary::getInstance()->writeLines(filenameVec, nodesVec, cellsVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeInterfaceCellsDebugCFC(Parameter *para)
{

    std::vector<UbTupleFloat3> nodesVec;
    std::vector<UbTupleInt8> cellsVec;
    int nodeNumberVec = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->coarseToFine.numberOfCells;
    }
    nodesVec.resize(nodeNumberVec * 8);
    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
            int pos = para->getParH(level)->coarseToFine.coarseCellIndices[u];

            double x1             = para->getParH(level)->coordinateX[pos];
            double x2             = para->getParH(level)->coordinateY[pos];
            double x3             = para->getParH(level)->coordinateZ[pos];
            double x1P            = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
            double x2P            = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
            double x3P            = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2P), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2P), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3P)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2), (float)(x3P)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2P), (float)(x3P)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2P), (float)(x3P)));

            cellsVec.push_back(makeUbTuple(nodeCount - 8, nodeCount - 7, nodeCount - 6, nodeCount - 5, nodeCount - 4,
                                           nodeCount - 3, nodeCount - 2, nodeCount - 1));
        }
        std::string filenameVec = para->getFName() + "_CellsCFC_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeOcts(filenameVec, nodesVec, cellsVec);
    }
}

//////////////////////////////////////////////////////////////////////////

void writeInterfaceCellsDebugCFF(Parameter *para)
{

    std::vector<UbTupleFloat3> nodesVec;
    std::vector<UbTupleInt8> cellsVec;
    int nodeNumberVec = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) // evtl. Maxlevel + 1
    {
        nodeNumberVec += (int)para->getParH(level)->coarseToFine.numberOfCells;
    }
    nodesVec.resize(nodeNumberVec * 8);
    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
            int pos = para->getParH(level)->coarseToFine.fineCellIndices[u];

            double x1             = para->getParH(level + 1)->coordinateX[pos];
            double x2             = para->getParH(level + 1)->coordinateY[pos];
            double x3             = para->getParH(level + 1)->coordinateZ[pos];
            double x1P            = para->getParH(level + 1)->coordinateX[para->getParH(level + 1)->neighborX[pos]];
            double x2P            = para->getParH(level + 1)->coordinateY[para->getParH(level + 1)->neighborY[pos]];
            double x3P            = para->getParH(level + 1)->coordinateZ[para->getParH(level + 1)->neighborZ[pos]];
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2P), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2P), (float)(x3)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3P)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2), (float)(x3P)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1P), (float)(x2P), (float)(x3P)));
            nodesVec[nodeCount++] = (makeUbTuple((float)(x1), (float)(x2P), (float)(x3P)));

            cellsVec.push_back(makeUbTuple(nodeCount - 8, nodeCount - 7, nodeCount - 6, nodeCount - 5, nodeCount - 4,
                                           nodeCount - 3, nodeCount - 2, nodeCount - 1));
        }
        std::string filenameVec = para->getFName() + "_CellsCFF_" + StringUtil::toString<int>(level);
        WbWriterVtkXmlBinary::getInstance()->writeOcts(filenameVec, nodesVec, cellsVec);
    }
}












//////////////////////////////////////////////////////////////////////////
// Functions for version with streams
//////////////////////////////////////////////////////////////////////////
void checkForSendOrRecvNode(int pos, int &commDir, int &commDirectionInCommAfterFtoC, int& indexInCommVector,
                            std::vector<ProcessNeighbor27> &sendRecvProcessNeighbor,
                            std::vector<ProcessNeighbor27> &sendRecvProcessNeighborsAfterFtoC, double indicator)
{
    for (uint pn = 0; pn < (uint)sendRecvProcessNeighbor.size(); pn++) {
        for (int j = 0; j < sendRecvProcessNeighbor[pn].numberOfNodes; j++) {
            if (pos == sendRecvProcessNeighbor[pn].index[j]) {
                commDir = indicator;
                indexInCommVector = j;
                if (j < sendRecvProcessNeighborsAfterFtoC[pn].numberOfNodes) {
                    commDirectionInCommAfterFtoC = indicator;
                }
                return;
            }
        }
    }
}

void checkForRecvNodeX(int pos, int &recvDir, int &recvDirectionInCommAfterFtoC, int& recvIndex, Parameter *para, int level)
{
    checkForSendOrRecvNode(pos, recvDir, recvDirectionInCommAfterFtoC, recvIndex, para->getParH(level)->recvProcessNeighborX,
                           para->getParH(level)->recvProcessNeighborsAfterFtoCX, 2.0);
}

void checkForRecvNodeY(int pos, int &recvDir, int &recvDirectionInCommAfterFtoC, int& recvIndex, Parameter *para, int level)
{
    checkForSendOrRecvNode(pos, recvDir, recvDirectionInCommAfterFtoC, recvIndex, para->getParH(level)->recvProcessNeighborY,
                           para->getParH(level)->recvProcessNeighborsAfterFtoCY, 4.0);
}

void checkForRecvNodeZ(int pos, int &recvDir, int &recvDirectionInCommAfterFtoC, int& recvIndex, Parameter *para, int level)
{
    checkForSendOrRecvNode(pos, recvDir, recvDirectionInCommAfterFtoC, recvIndex, para->getParH(level)->recvProcessNeighborZ,
                           para->getParH(level)->recvProcessNeighborsAfterFtoCZ, 8.0);
}

void checkForSendNodeX(int pos, int &sendDir, int &sendDirectionInCommAfterFtoC, int& sendIndex, Parameter *para, int level)
{
    checkForSendOrRecvNode(pos, sendDir, sendDirectionInCommAfterFtoC, sendIndex, para->getParH(level)->sendProcessNeighborX,
                           para->getParH(level)->sendProcessNeighborsAfterFtoCX, 2.0);
}

void checkForSendNodeY(int pos, int &sendDir, int &sendDirectionInCommAfterFtoC, int& sendIndex, Parameter *para, int level)
{
    checkForSendOrRecvNode(pos, sendDir, sendDirectionInCommAfterFtoC, sendIndex, para->getParH(level)->sendProcessNeighborY,
                           para->getParH(level)->sendProcessNeighborsAfterFtoCY, 4.0);
}

void checkForSendNodeZ(int pos, int &sendDir, int &sendDirectionInCommAfterFtoC, int& sendIndex, Parameter *para, int level)
{
    checkForSendOrRecvNode(pos, sendDir, sendDirectionInCommAfterFtoC, sendIndex, para->getParH(level)->sendProcessNeighborZ,
                           para->getParH(level)->sendProcessNeighborsAfterFtoCZ, 8.0);
}

void writeInterfaceFCC_Send(Parameter *para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    // nodedata
    std::vector<std::string> datanames = { "sparse index", "borderBulk", "sendDirection",
                                           "sendDirectionInCommAfterFtoC", "sendIndex" };
    // sendDirection: x = 2, y = 4, z = 8
    // borderBulk: border = 1, bulk = 0
    std::vector<std::vector<double>> nodedata;

    for (int level = 0; level < para->getMaxLevel(); level++) {
        nodeNumberVec += (int)para->getParH(level)->fineToCoarse.numberOfCells;
    }

    nodesVec.resize(nodeNumberVec);
    nodedata.resize(datanames.size(), std::vector<double>(nodeNumberVec));

    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++) {
            int pos                = para->getParH(level)->fineToCoarse.coarseCellIndices[u];
            nodedata[0][nodeCount] = pos;

            // coordinate section
            double x1           = para->getParH(level)->coordinateX[pos];
            double x2           = para->getParH(level)->coordinateY[pos];
            double x3           = para->getParH(level)->coordinateZ[pos];
            nodesVec[nodeCount] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));

            // nodedata section
            nodedata[1][nodeCount]           = u < para->getParH(level)->fineToCoarseBorder.numberOfCells;
            int sendDir                      = 0.0;
            int sendDirectionInCommAfterFtoC = 0.0;
            int sendIndex                    = 0.0;

            checkForSendNodeX(pos, sendDir, sendIndex, sendDirectionInCommAfterFtoC, para, level);
            checkForSendNodeY(pos, sendDir, sendIndex, sendDirectionInCommAfterFtoC, para, level);
            checkForSendNodeZ(pos, sendDir, sendIndex, sendDirectionInCommAfterFtoC, para, level);
            nodedata[2][nodeCount] = sendDir;
            nodedata[3][nodeCount] = sendDirectionInCommAfterFtoC;
            nodedata[4][nodeCount] = sendIndex;

            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeInterfaceFCC_Send_PID_" +
                                  std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}

void writeInterfaceCFC_Recv(Parameter *para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;
    int nodeNumberVec = 0;

    // nodedata
    std::vector<std::string> datanames = { "sparse index", "borderBulk", "recvDirection",
                                           "recvDirectionInCommAfterFtoC", "recvIndex"};
    // recvDirection: x = 2, y = 4, z = 8
    // borderBulk: border = 1, bulk = 0
    std::vector<std::vector<double>> nodedata;

    for (int level = 0; level < para->getMaxLevel(); level++) {
        nodeNumberVec += (int)para->getParH(level)->coarseToFine.numberOfCells;
    }

    nodesVec.resize(nodeNumberVec);
    nodedata.resize(datanames.size(), std::vector<double>(nodeNumberVec));

    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++) {
            int pos                = para->getParH(level)->coarseToFine.coarseCellIndices[u];
            nodedata[0][nodeCount] = pos;

            // coordinate section
            double x1           = para->getParH(level)->coordinateX[pos];
            double x2           = para->getParH(level)->coordinateY[pos];
            double x3           = para->getParH(level)->coordinateZ[pos];
            nodesVec[nodeCount] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));

            // nodedata section
            nodedata[1][nodeCount]           = u < para->getParH(level)->coarseToFineBorder.numberOfCells;
            int recvDir                      = 0.0;
            int recvDirectionInCommAfterFtoC = 0.0;
            int recvIndex                    = 0.0;

            checkForRecvNodeX(pos, recvDir, recvIndex, recvDirectionInCommAfterFtoC, para, level);
            checkForRecvNodeY(pos, recvDir, recvIndex, recvDirectionInCommAfterFtoC, para, level);
            checkForRecvNodeZ(pos, recvDir, recvIndex, recvDirectionInCommAfterFtoC, para, level);
            nodedata[2][nodeCount] = recvDir;
            nodedata[3][nodeCount] = recvDirectionInCommAfterFtoC;
            nodedata[4][nodeCount] = recvIndex;
            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeInterfaceCFC_Recv_PID_" +
                                  std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}

void addToNodesVector(const int level, const int pos, std::vector<UbTupleFloat3> &nodesVec, Parameter *para)
{
    double x1 = para->getParH(level)->coordinateX[pos];
    double x2 = para->getParH(level)->coordinateY[pos];
    double x3 = para->getParH(level)->coordinateZ[pos];
    nodesVec.push_back(makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
}

void writeSendNodesStream(Parameter *para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;

    // nodedata
    std::vector<std::string> datanames = { "sparse index", "sendDirection", "sendDirectionInCommAfterFtoC", "sendIndex",
                                           "inICcellFCC" };
    // sendDirection: x = 2, y = 4, z = 8
    std::vector<std::vector<double>> nodedata;
    nodedata.resize(datanames.size());

    int pos;
    int sendDirectionInCommAfterFtoC;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        // X
        for (int pn = 0; pn < (int)para->getParH(level)->sendProcessNeighborX.size(); pn++) {
            for (int i = 0; i < para->getParH(level)->sendProcessNeighborX[pn].numberOfNodes; i++) {
                pos = para->getParH(level)->sendProcessNeighborX[pn].index[i];
                nodedata[0].push_back(pos);
                addToNodesVector(level, pos, nodesVec, para);

                nodedata[1].push_back(2.0);
                sendDirectionInCommAfterFtoC =
                    (i < para->getParH(level)->sendProcessNeighborsAfterFtoCX[pn].numberOfNodes) ? 2.0 : 0.0;
                nodedata[2].push_back(sendDirectionInCommAfterFtoC);
                nodedata[3].push_back((double)i);
            }
        }

        // Y
        for (int pn = 0; pn < (int)para->getParH(level)->sendProcessNeighborY.size(); pn++) {
            for (int i = 0; i < para->getParH(level)->sendProcessNeighborY[pn].numberOfNodes; i++) {
                pos = para->getParH(level)->sendProcessNeighborY[pn].index[i];

                sendDirectionInCommAfterFtoC =
                    (i < para->getParH(level)->sendProcessNeighborsAfterFtoCY[pn].numberOfNodes) ? 4.0 : 0.0;

                auto it = std::find(nodedata[0].begin(), nodedata[0].end(), pos);
                if (it == nodedata[0].end()) {
                    nodedata[0].push_back(pos);
                    addToNodesVector(level, pos, nodesVec, para);
                    nodedata[1].push_back(4.0);
                    nodedata[2].push_back(sendDirectionInCommAfterFtoC);
                    nodedata[3].push_back((double) i);
                } else {
                    int posInVectors = it - nodedata[0].begin();
                    nodedata[1][posInVectors] += 4.0;
                    nodedata[2][posInVectors] += sendDirectionInCommAfterFtoC;
                    nodedata[3][posInVectors] = (double)i;
                }
            }
        }

        // Z
        for (int pn = 0; pn < (int)para->getParH(level)->sendProcessNeighborZ.size(); pn++) {
            for (int i = 0; i < para->getParH(level)->sendProcessNeighborZ[pn].numberOfNodes; i++) {
                pos = para->getParH(level)->sendProcessNeighborZ[pn].index[i];

                sendDirectionInCommAfterFtoC =
                    (i < para->getParH(level)->sendProcessNeighborsAfterFtoCZ[pn].numberOfNodes) ? 8.0 : 0.0;

                auto it = std::find(nodedata[0].begin(), nodedata[0].end(), pos);
                if (it == nodedata[0].end()) {
                    nodedata[0].push_back(pos);
                    addToNodesVector(level, pos, nodesVec, para);
                    nodedata[1].push_back(8.0);
                    nodedata[2].push_back(sendDirectionInCommAfterFtoC);
                    nodedata[3].push_back((double) i);
                } else {
                    int posInVectors = it - nodedata[0].begin();
                    nodedata[1][posInVectors] += 8.0;
                    nodedata[2][posInVectors] += sendDirectionInCommAfterFtoC;
                    nodedata[3][posInVectors] = (double)i;
                }
            }
        }

        // check if node is in a coarse cell for the interpolation from fine to coarse
        nodedata[4].resize(nodedata[0].size());
        for (int i = 0; i < (int)nodedata[0].size(); i++) {
            pos = nodedata[0][i];
            for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++) {
                if (para->getParH(level)->fineToCoarse.coarseCellIndices[u] == (uint)pos) {
                    nodedata[4][i] = 1.0;
                    break;
                }
                nodedata[4][i] = 0.0;
            }
        }
        std::string filenameVec = para->getFName() + "_writeSendNodesStreams_PID_" +
                                  std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}

void writeRecvNodesStream(Parameter *para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;

    // nodedata
    std::vector<std::string> datanames = { "sparse index", "recvDirection", "recvDirectionInCommAfterFtoC", "recvIndex" };
    // sendDirection: x = 2, y = 4, z = 8
    std::vector<std::vector<double>> nodedata;
    nodedata.resize(datanames.size());

    int pos;
    int recvDirectionInCommAfterFtoC;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        // X
        for (int pn = 0; pn < (int)para->getParH(level)->recvProcessNeighborX.size(); pn++) {
            for (int i = 0; i < para->getParH(level)->recvProcessNeighborX[pn].numberOfNodes; i++) {
                pos = para->getParH(level)->recvProcessNeighborX[pn].index[i];
                nodedata[0].push_back(pos);
                addToNodesVector(level, pos, nodesVec, para);

                nodedata[1].push_back(2.0);
                recvDirectionInCommAfterFtoC =
                    (i < para->getParH(level)->recvProcessNeighborsAfterFtoCX[pn].numberOfNodes) ? 2.0 : 0.0;
                nodedata[2].push_back(recvDirectionInCommAfterFtoC);
                nodedata[3].push_back(i);
            }
        }

        // Y
        for (int pn = 0; pn < (int)para->getParH(level)->recvProcessNeighborY.size(); pn++) {
            for (int i = 0; i < para->getParH(level)->recvProcessNeighborY[pn].numberOfNodes; i++) {
                pos = para->getParH(level)->recvProcessNeighborY[pn].index[i];

                recvDirectionInCommAfterFtoC =
                    (i < para->getParH(level)->recvProcessNeighborsAfterFtoCY[pn].numberOfNodes) ? 4.0 : 0.0;

                auto it = std::find(nodedata[0].begin(), nodedata[0].end(), pos);
                if (it == nodedata[0].end()) {
                    nodedata[0].push_back(pos);
                    addToNodesVector(level, pos, nodesVec, para);
                    nodedata[1].push_back(4.0);
                    nodedata[2].push_back(recvDirectionInCommAfterFtoC);
                    nodedata[3].push_back(i);
                } else {
                    int posInVectors = it - nodedata[0].begin();
                    nodedata[1][posInVectors] += 4.0;
                    nodedata[2][posInVectors] += recvDirectionInCommAfterFtoC;
                    nodedata[3][posInVectors] += i;
                }
            }
        }

        // Z
        for (int pn = 0; pn < (int)para->getParH(level)->recvProcessNeighborZ.size(); pn++) {
            for (int i = 0; i < para->getParH(level)->recvProcessNeighborZ[pn].numberOfNodes; i++) {
                pos = para->getParH(level)->recvProcessNeighborZ[pn].index[i];

                recvDirectionInCommAfterFtoC =
                    (i < para->getParH(level)->recvProcessNeighborsAfterFtoCZ[pn].numberOfNodes) ? 8.0 : 0.0;

                auto it = std::find(nodedata[0].begin(), nodedata[0].end(), pos);
                if (it == nodedata[0].end()) {
                    nodedata[0].push_back(pos);
                    addToNodesVector(level, pos, nodesVec, para);
                    nodedata[1].push_back(8.0);
                    nodedata[2].push_back(recvDirectionInCommAfterFtoC);
                    nodedata[3].push_back(i);
                } else {
                    int posInVectors = it - nodedata[0].begin();
                    nodedata[1][posInVectors] += 8.0;
                    nodedata[2][posInVectors] += recvDirectionInCommAfterFtoC;
                    nodedata[3][posInVectors] += i;
                }
            }
        }

        // Recv are nodes ghost nodes and therefore they can't be coarse cells for the interpolation from coarse to fine

        std::string filenameVec = para->getFName() + "_writeRecvNodesStreams_PID_" +
                                  std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}

} // namespace InterfaceDebugWriter
#endif
