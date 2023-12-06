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
#ifndef UNSTRUCTUREDGRID_HPP
#define UNSTRUCTUREDGRID_HPP

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <lbm/constants/D3Q27.h>

#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

namespace UnstructuredGridWriter
{

void writeUnstructuredGrid(Parameter* para, int level, std::string& fname, std::string& filenameVec2)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleUInt8 > cells;
    std::vector< string > nodedatanames;
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    //int number1,number2,number3,number4,number5,number6,number7,number8;
    vector< vector< double > > nodedata(nodedatanames.size());

    bool neighborsFluid;

    unsigned long long allnodes = para->getParH(level)->numberOfNodes * 8;

    nodes.resize(allnodes);
    nodedata[0].resize(allnodes);
    nodedata[1].resize(allnodes);
    nodedata[2].resize(allnodes);
    nodedata[3].resize(allnodes);
    nodedata[4].resize(allnodes);

    unsigned int nodeCount = 0;
    double nodeDeltaLevel = para->getParH(level)->dx;

    for (size_t pos = 0; pos < para->getParH(level)->numberOfNodes; pos++)
    {
        if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID /*!= GEO_VOID*/)
        {
            //////////////////////////////////////////////////////////////////////////
            double ix1  = para->getParH(level)->coordinateX[pos];//-STARTOFFX;
            double ix2  = para->getParH(level)->coordinateY[pos];//-STARTOFFY;
            double ix3  = para->getParH(level)->coordinateZ[pos];//-STARTOFFZ;
            double ix1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];//-STARTOFFX;
            double ix2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];//-STARTOFFY;
            double ix3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];//-STARTOFFZ;
            //////////////////////////////////////////////////////////////////////////
            double x1  = ix1;  // para->getParH(level)->distX + ix1 *nodeDeltaLevel;// + tmpDist;
            double x2  = ix2;  // para->getParH(level)->distY + ix2 *nodeDeltaLevel;// + tmpDist;
            double x3  = ix3;  // para->getParH(level)->distZ + ix3 *nodeDeltaLevel;// + tmpDist;
            double x1P = ix1P; // para->getParH(level)->distX + ix1P*nodeDeltaLevel;// + tmpDist;
            double x2P = ix2P; // para->getParH(level)->distY + ix2P*nodeDeltaLevel;// + tmpDist;
            double x3P = ix3P; // para->getParH(level)->distZ + ix3P*nodeDeltaLevel;// + tmpDist;
            //////////////////////////////////////////////////////////////////////////
            neighborsFluid = true;
            //////////////////////////////////////////////////////////////////////////
            //1
            nodes[nodeCount]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[pos] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[pos] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[pos] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[pos];
            //if(para->getParH(level)->typeOfGridNode[pos]==GEO_VOID) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //2
            nodes[nodeCount]=( makeUbTuple( (float)(x1P),(float)(x2 ),(float)(x3 ) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborX[pos]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborX[pos]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborX[pos]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborX[pos]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborX[pos]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborX[pos]]==GEO_VOID) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //3
            nodes[nodeCount]=( makeUbTuple( (float)(x1P),(float)(x2P),(float)(x3 ) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]==GEO_VOID) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //4
            nodes[nodeCount]=( makeUbTuple( (float)(x1 ),(float)(x2P),(float)(x3 ) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborY[pos]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborY[pos]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborY[pos]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborY[pos]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborY[pos]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborY[pos]]==GEO_VOID) neighborsFluid = false;
            //if((para->getParH(level)->neighborY[pos]<=pos) && ((para->getParH(level)->coordinateY[pos]) > (para->getParH(level)->gridNY-2))) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //5
            nodes[nodeCount]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3P) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborZ[pos]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborZ[pos]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborZ[pos]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborZ[pos]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[pos]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[pos]]==GEO_VOID) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //6
            nodes[nodeCount]=( makeUbTuple( (float)(x1P),(float)(x2 ),(float)(x3P) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborZ[para->getParH(level)->neighborX[pos]]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborZ[para->getParH(level)->neighborX[pos]]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborZ[para->getParH(level)->neighborX[pos]]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborZ[para->getParH(level)->neighborX[pos]]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[para->getParH(level)->neighborX[pos]]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[para->getParH(level)->neighborX[pos]]]==GEO_VOID) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //7
            nodes[nodeCount]=( makeUbTuple( (float)(x1P),(float)(x2P),(float)(x3P) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[para->getParH(level)->neighborX[pos]]]]==GEO_VOID) neighborsFluid = false;
            nodeCount++;
            //////////////////////////////////////////////////////////////////////////
            //8
            nodes[nodeCount]=( makeUbTuple( (float)(x1 ),(float)(x2P),(float)(x3P) ) );
            nodedata[0][nodeCount] = para->getParH(level)->rho[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[pos]]] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][nodeCount] = para->getParH(level)->velocityX[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[pos]]] * para->getVelocityRatio();
            nodedata[2][nodeCount] = para->getParH(level)->velocityY[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[pos]]] * para->getVelocityRatio();
            nodedata[3][nodeCount] = para->getParH(level)->velocityZ[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[pos]]] * para->getVelocityRatio();
            nodedata[4][nodeCount] = para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[pos]]];
            //if(para->getParH(level)->typeOfGridNode[para->getParH(level)->neighborZ[para->getParH(level)->neighborY[pos]]]==GEO_VOID) neighborsFluid = false;
            nodeCount++;

            if(neighborsFluid)
            {
                cells.push_back( makeUbTuple(nodeCount-8,nodeCount-7,nodeCount-6,nodeCount-5,nodeCount-4,nodeCount-3,nodeCount-2,nodeCount-1) );        
            }
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
    //WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec2,nodes);
}
//////////////////////////////////////////////////////////////////////////



bool isPeriodicCell(Parameter* para, int level, unsigned int number2, unsigned int number1, unsigned int number3, unsigned int number5)
{
    return (para->getParH(level)->coordinateX[number2] < para->getParH(level)->coordinateX[number1]) ||
        (para->getParH(level)->coordinateY[number3] < para->getParH(level)->coordinateY[number1]) ||
        (para->getParH(level)->coordinateZ[number5] < para->getParH(level)->coordinateZ[number1]);
}


//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridLT(Parameter* para, int level, vector<string >& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());


    //printf("\n test for if... \n");
    for (unsigned int part=0; part < fname.size(); part++)
    {
        vxmax = 0;
        // printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////

        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);

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
        //int counter = 0;
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");

        for (unsigned int pos=startpos;pos<endpos;pos++)
        {
            if (/*para->getParH(level)->typeOfGridNode[pos] >= GEO_FLUID*/true)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor node data... \n");
                nodes[dn1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor numbers... \n");
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor neighborsFluid... \n");
                if (para->getParH(level)->typeOfGridNode[number2] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] < GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                //if(neighborsFluid==false) counter++;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor numbers and neighborsFluid... \n");
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                //if(neighborsFluid==false) counter++;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor dn... \n");
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                //if( std::fabs(nodedata[2][dn1]) > std::fabs(vxmax) ) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                //counter++;
                if (neighborsFluid) cells.push_back( makeUbTuple(dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8) );
                //////////////////////////////////////////////////////////////////////////
            }
            //printf("\n test II... \n");
        }
        //printf("\n number of cells: %d at level: %d\n", cells.size(), level);
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part],nodes,cells,nodedatanames,nodedata);
        //WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(fname[part], nodes, nodedatanames, nodedata);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n vx max: %.1f at level: %d\n", vxmax, level);
        //printf("\n counter: %d at level: %d\n", counter, level);
    } 
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridLTwithTurbulentViscosity(Parameter* para, int level, vector<string >& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    nodedatanames.push_back("turbVis");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());


    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////

        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);

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
        nodedata[6].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //int counter = 0;
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");

        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (/*para->getParH(level)->typeOfGridNode[pos] >= GEO_FLUID*/true)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                nodedata[6][dn1] = (double)para->getParH(level)->turbViscosity[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] < GEO_FLUID)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
        //WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(fname[part], nodes, nodedatanames, nodedata);
    }
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridLTwithTurbulentViscosityDebug(Parameter* para, int level, vector<string >& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    nodedatanames.push_back("turbVis");
    nodedatanames.push_back("gSij ");
    nodedatanames.push_back("gSDij");
    nodedatanames.push_back("gDxvx");
    nodedatanames.push_back("gDyvx");
    nodedatanames.push_back("gDzvx");
    nodedatanames.push_back("gDxvy");
    nodedatanames.push_back("gDyvy");
    nodedatanames.push_back("gDzvy");
    nodedatanames.push_back("gDxvz");
    nodedatanames.push_back("gDyvz");
    nodedatanames.push_back("gDzvz");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());


    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        nodedata[6].resize(sizeOfNodes);
        nodedata[7].resize(sizeOfNodes);
        nodedata[8].resize(sizeOfNodes);
        nodedata[9].resize(sizeOfNodes);
        nodedata[10].resize(sizeOfNodes);
        nodedata[11].resize(sizeOfNodes);
        nodedata[12].resize(sizeOfNodes);
        nodedata[13].resize(sizeOfNodes);
        nodedata[14].resize(sizeOfNodes);
        nodedata[15].resize(sizeOfNodes);
        nodedata[16].resize(sizeOfNodes);
        nodedata[17].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //int counter = 0;
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");

        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (/*para->getParH(level)->typeOfGridNode[pos] >= GEO_FLUID*/true)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                nodedata[6][dn1] = (double)para->getParH(level)->turbViscosity[pos] * (double)para->getViscosityRatio();
                nodedata[7][dn1] = (double)para->getParH(level)->gSij[pos] * (double)para->getVelocityRatio();
                nodedata[8][dn1] = (double)para->getParH(level)->gSDij[pos] * (double)para->getVelocityRatio();
                nodedata[9][dn1] = (double)para->getParH(level)->gDxvx[pos] * (double)para->getVelocityRatio();
                nodedata[10][dn1] = (double)para->getParH(level)->gDyvx[pos] * (double)para->getVelocityRatio();
                nodedata[11][dn1] = (double)para->getParH(level)->gDzvx[pos] * (double)para->getVelocityRatio();
                nodedata[12][dn1] = (double)para->getParH(level)->gDxvy[pos] * (double)para->getVelocityRatio();
                nodedata[13][dn1] = (double)para->getParH(level)->gDyvy[pos] * (double)para->getVelocityRatio();
                nodedata[14][dn1] = (double)para->getParH(level)->gDzvy[pos] * (double)para->getVelocityRatio();
                nodedata[15][dn1] = (double)para->getParH(level)->gDxvz[pos] * (double)para->getVelocityRatio();
                nodedata[16][dn1] = (double)para->getParH(level)->gDyvz[pos] * (double)para->getVelocityRatio();
                nodedata[17][dn1] = (double)para->getParH(level)->gDzvz[pos] * (double)para->getVelocityRatio();
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] < GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] < GEO_FLUID)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
        //WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(fname[part], nodes, nodedatanames, nodedata);
    }
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridPM(Parameter* para, int level, vector<string >& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());


    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        //int counter = 0;
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");

        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if ((para->getParH(level)->typeOfGridNode[pos] >= GEO_FLUID) || ((para->getParH(level)->typeOfGridNode[pos] >= GEO_PM_0) && (para->getParH(level)->typeOfGridNode[pos] <= GEO_PM_2)))
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor node data... \n");
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor numbers... \n");
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor neighborsFluid... \n");
                if (((para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number2] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number2] > GEO_PM_2)) ||
                    ((para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number3] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number3] > GEO_PM_2)) ||
                    ((para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number4] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number4] > GEO_PM_2)) ||
                    ((para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number5] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number5] > GEO_PM_2)) ||
                    ((para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number6] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number6] > GEO_PM_2)) ||
                    ((para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number7] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number7] > GEO_PM_2)) ||
                    ((para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID) && (para->getParH(level)->typeOfGridNode[number8] < GEO_PM_0) && (para->getParH(level)->typeOfGridNode[number8] > GEO_PM_2)))  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                //if(neighborsFluid==false) counter++;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor numbers and neighborsFluid... \n");
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                //if(neighborsFluid==false) counter++;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test vor dn... \n");
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                //if( std::fabs(nodedata[2][dn1]) > std::fabs(vxmax) ) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////

                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;

                //counter++;
                if (neighborsFluid) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
            //printf("\n test II... \n");
        }
        //printf("\n number of cells: %d at level: %d\n", cells.size(), level);
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
        //WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(fname[part], nodes, nodedatanames, nodedata);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n vx max: %.1f at level: %d\n", vxmax, level);
        //printf("\n counter: %d at level: %d\n", counter, level);
    }
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridLTConc(Parameter* para, int level, vector<string >& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    nodedatanames.push_back("Conc");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part=0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        nodedata[6].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");
        for (unsigned int pos=startpos;pos<endpos;pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                nodedata[6][dn1] = (double)para->getParH(level)->concentration[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                //if( std::fabs(nodedata[2][dn1]) > std::fabs(vxmax) ) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid) cells.push_back( makeUbTuple(dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part],nodes,cells,nodedatanames,nodedata);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n vx max: %.1f at level: %d\n", vxmax, level);
    } 
}
//////////////////////////////////////////////////////////////////////////











//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridBig(Parameter* para, int level, std::string& fname, std::string& fname2) 
{
    unsigned int limitOfNodes = 30000000; //27 Million
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    //double posmax = 0;
    double vxmax = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    if ((uint)para->getParH(level)->numberOfNodes > limitOfNodes)
    {
        //printf("\n test in if I... \n");
        unsigned int restOfNodes = (uint)para->getParH(level)->numberOfNodes - limitOfNodes;
        //////////////////////////////////////////////////////////////////////////
        //PART I
        nodes.resize(limitOfNodes);
        nodedata[0].resize(limitOfNodes);
        nodedata[1].resize(limitOfNodes);
        nodedata[2].resize(limitOfNodes);
        nodedata[3].resize(limitOfNodes);
        nodedata[4].resize(limitOfNodes);
        nodedata[5].resize(limitOfNodes);

        //printf("\n test in if II... \n");
        for (unsigned int pos=0;pos<limitOfNodes;pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][number1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][number1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][number1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][number1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][number1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][number1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > limitOfNodes ||
                    number3 > limitOfNodes ||
                    number4 > limitOfNodes ||
                    number5 > limitOfNodes ||
                    number6 > limitOfNodes ||
                    number7 > limitOfNodes ||
                    number8 > limitOfNodes )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                //if (level == 0 &&
                //    (number2 <= number1 ||
                //    number3 <= number1 ||
                //    number4 <= number1 ||
                //    number5 <= number1 ||
                //    number6 <= number1 ||
                //    number7 <= number1 ||
                //    number8 <= number1) )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if( std::fabs(nodedata[2][number1]) > std::fabs(vxmax) ) vxmax = nodedata[2][number1];
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);

        //printf("\n test in if III... \n");
        //////////////////////////////////////////////////////////////////////////
        //PART II
        nodes.resize(restOfNodes);
        nodedata[0].resize(restOfNodes);
        nodedata[1].resize(restOfNodes);
        nodedata[2].resize(restOfNodes);
        nodedata[3].resize(restOfNodes);
        nodedata[4].resize(restOfNodes);
        nodedata[5].resize(restOfNodes);
        //printf("\n test in if IV... \n");

        for (size_t pos = limitOfNodes; pos < para->getParH(level)->numberOfNodes; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - limitOfNodes;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - limitOfNodes;
                dn3 = number3 - limitOfNodes;
                dn4 = number4 - limitOfNodes;
                dn5 = number5 - limitOfNodes;
                dn6 = number6 - limitOfNodes;
                dn7 = number7 - limitOfNodes;
                dn8 = number8 - limitOfNodes;
                //////////////////////////////////////////////////////////////////////////
                if( std::fabs(nodedata[2][dn1]) > std::fabs(vxmax) ) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid==true) cells2.push_back( makeUbTuple(dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        //printf("\n test in if V... \n");
        //WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeDataMS(fname,nodes,cells2,nodedatanames,nodedata);
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname2,nodes,cells2,nodedatanames,nodedata);
        //printf("\n test in if VI... \n");
        //////////////////////////////////////////////////////////////////////////
        //printf("pos max: %.1f", posmax);
        printf("\n vx max: %.1f at level: %d\n", vxmax, level);
    } 
    else
    {
        //printf("\n test in else I... \n");
        nodes.resize(para->getParH(level)->numberOfNodes);
        nodedata[0].resize(para->getParH(level)->numberOfNodes);
        nodedata[1].resize(para->getParH(level)->numberOfNodes);
        nodedata[2].resize(para->getParH(level)->numberOfNodes);
        nodedata[3].resize(para->getParH(level)->numberOfNodes);
        nodedata[4].resize(para->getParH(level)->numberOfNodes);
        nodedata[5].resize(para->getParH(level)->numberOfNodes);

        //printf("\n test in else II... \n");
        for (size_t pos = 0; pos < para->getParH(level)->numberOfNodes; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //printf("\n test in else-for I pos = %d \n", pos);
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test in else-for II pos = %d \n", pos);
                nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][number1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][number1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][number1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][number1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][number1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][number1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test in else-for III pos = %d \n", pos);
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test in else-for VI pos = %d \n", pos);
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                //if (level == 0 &&
                //    (number2 <= number1 ||
                //    number3 <= number1 ||
                //    number4 <= number1 ||
                //    number5 <= number1 ||
                //    number6 <= number1 ||
                //    number7 <= number1 ||
                //    number8 <= number1) )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if( std::fabs(nodedata[2][number1]) > std::fabs(vxmax) ) vxmax = nodedata[2][number1];
                //////////////////////////////////////////////////////////////////////////
                //printf("\n test in else-for V pos = %d \n", pos);
                if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        //printf("\n test in else III... \n");
        //WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeDataMS(fname,nodes,cells,nodedatanames,nodedata);
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
        //WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec2,nodes);
        //printf("\n test in else IV... \n");
        //////////////////////////////////////////////////////////////////////////
        //printf("pos max: %.1f", posmax);
        printf("\n vx max: %.1f at level: %d\n", vxmax, level);
    }
}
//////////////////////////////////////////////////////////////////////////











//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridEff(Parameter* para, int level, std::string& fname, std::string& filenameVec2) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    bool neighborsFluid;
    double vxmax = 0;
    vector<vector<double>> nodedata(nodedatanames.size());

    nodes.resize(para->getParH(level)->numberOfNodes);
    nodedata[0].resize(para->getParH(level)->numberOfNodes);
    nodedata[1].resize(para->getParH(level)->numberOfNodes);
    nodedata[2].resize(para->getParH(level)->numberOfNodes);
    nodedata[3].resize(para->getParH(level)->numberOfNodes);
    nodedata[4].resize(para->getParH(level)->numberOfNodes);
    nodedata[5].resize(para->getParH(level)->numberOfNodes);

    for (size_t pos = 0; pos < para->getParH(level)->numberOfNodes; pos++)
    {
        if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
        {
            //////////////////////////////////////////////////////////////////////////
            double x1  = para->getParH(level)->coordinateX[pos];
            double x2  = para->getParH(level)->coordinateY[pos];
            double x3  = para->getParH(level)->coordinateZ[pos];
            double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
            double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
            double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
            //////////////////////////////////////////////////////////////////////////
            number1 = pos;
            neighborsFluid = true;
            //////////////////////////////////////////////////////////////////////////
            nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
            nodedata[0][number1] = (double)pos;
            //nodedata[0][number1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
            nodedata[1][number1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
            nodedata[2][number1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
            nodedata[3][number1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
            nodedata[4][number1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
            nodedata[5][number1] = (double)para->getParH(level)->typeOfGridNode[pos];
            //////////////////////////////////////////////////////////////////////////
            number2 = para->getParH(level)->neighborX[number1];
            number3 = para->getParH(level)->neighborY[number2];
            number4 = para->getParH(level)->neighborY[number1];
            number5 = para->getParH(level)->neighborZ[number1];
            number6 = para->getParH(level)->neighborZ[number2];
            number7 = para->getParH(level)->neighborZ[number3];
            number8 = para->getParH(level)->neighborZ[number4];
            //////////////////////////////////////////////////////////////////////////
            if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            //if (level == 0 &&
            //    (number2 <= number1 ||
            //    number3 <= number1 ||
            //    number4 <= number1 ||
            //    number5 <= number1 ||
            //    number6 <= number1 ||
            //    number7 <= number1 ||
            //    number8 <= number1) )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
            //////////////////////////////////////////////////////////////////////////
        }
        if( std::fabs(nodedata[2][number1]) > std::fabs(vxmax) ) vxmax = nodedata[2][number1];
    }
    WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
    printf("\n vx max: %.1f at level: %d\n", vxmax, level);
}




//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridAsciiEff(Parameter* para, int level, std::string& fname, std::string& filenameVec2) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    vector< string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    bool neighborsFluid;
    double posmax = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    nodes.resize(para->getParH(level)->numberOfNodes);
    nodedata[0].resize(para->getParH(level)->numberOfNodes);
    nodedata[1].resize(para->getParH(level)->numberOfNodes);
    nodedata[2].resize(para->getParH(level)->numberOfNodes);
    nodedata[3].resize(para->getParH(level)->numberOfNodes);
    nodedata[4].resize(para->getParH(level)->numberOfNodes);
    nodedata[5].resize(para->getParH(level)->numberOfNodes);

    for (size_t pos = 0; pos < para->getParH(level)->numberOfNodes; pos++)
    {
        if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
        {
            //////////////////////////////////////////////////////////////////////////
            double x1  = para->getParH(level)->coordinateX[pos];
            double x2  = para->getParH(level)->coordinateY[pos];
            double x3  = para->getParH(level)->coordinateZ[pos];
            double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
            double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
            double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
            //////////////////////////////////////////////////////////////////////////
            number1 = pos;
            neighborsFluid = true;
            //////////////////////////////////////////////////////////////////////////
            nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
            nodedata[0][number1] = (double)pos;
            //nodedata[0][number1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
            nodedata[1][number1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
            nodedata[2][number1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
            nodedata[3][number1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
            nodedata[4][number1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
            nodedata[5][number1] = (double)para->getParH(level)->typeOfGridNode[pos];
            //////////////////////////////////////////////////////////////////////////
            number2 = para->getParH(level)->neighborX[number1];
            number3 = para->getParH(level)->neighborY[number2];
            number4 = para->getParH(level)->neighborY[number1];
            number5 = para->getParH(level)->neighborZ[number1];
            number6 = para->getParH(level)->neighborZ[number2];
            number7 = para->getParH(level)->neighborZ[number3];
            number8 = para->getParH(level)->neighborZ[number4];
            //////////////////////////////////////////////////////////////////////////
            if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            //if (level == 0 &&
            //    (number2 <= number1 ||
            //    number3 <= number1 ||
            //    number4 <= number1 ||
            //    number5 <= number1 ||
            //    number6 <= number1 ||
            //    number7 <= number1 ||
            //    number8 <= number1) )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
            //////////////////////////////////////////////////////////////////////////
        }
        if( std::fabs((double)pos) > std::fabs(posmax) ) posmax = (double)pos;
    }
    WbWriterVtkXmlASCII::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
    //WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
    ////WbWriterVtkXmlBinary::getInstance()->writeNodes(filenameVec2,nodes);
    printf("\ncells: %.1f \n", (double)cells.size());
    printf("nodes: %.1f \n", (double)nodes.size());
    printf("pos max: %.1f \n", posmax);
}
//////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridMeanLT(Parameter* para, int level, vector<string >& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("pressMed");
    nodedatanames.push_back("rhoMed");
    nodedatanames.push_back("vx1Med");
    nodedatanames.push_back("vx2Med");
    nodedatanames.push_back("vx3Med");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part=0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////

        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);

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
        //printf("\n test in if II... \n");
        for (size_t pos = startpos; pos < endpos; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][dn1] = para->getParH(level)->meanPressureOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[1][dn1] = para->getParH(level)->meanDensityOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[2][dn1] = para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio();
                nodedata[3][dn1] = para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio();
                nodedata[4][dn1] = para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if( std::fabs(nodedata[2][dn1]) > std::fabs(vxmax) ) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid==true) cells.push_back( makeUbTuple(dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part],nodes,cells,nodedatanames,nodedata);
        //////////////////////////////////////////////////////////////////////////
        printf("\n vx mean max: %.1f at level: %d\n", vxmax, level);
    } 
}
//////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridMeanLTConc(Parameter* para, int level, vector<string >& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("concMed");
    nodedatanames.push_back("pressMed");
    nodedatanames.push_back("rhoMed");
    nodedatanames.push_back("vx1Med");
    nodedatanames.push_back("vx2Med");
    nodedatanames.push_back("vx3Med");
    nodedatanames.push_back("geo");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        //printf("\n test in if II... \n");
        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = para->getParH(level)->meanConcentrationOut[pos];
                nodedata[1][dn1] = para->getParH(level)->meanPressureOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[2][dn1] = para->getParH(level)->meanDensityOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[3][dn1] = para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio();
                nodedata[4][dn1] = para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio();
                nodedata[5][dn1] = para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio();
                nodedata[6][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if (std::fabs(nodedata[2][dn1]) > std::fabs(vxmax)) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid == true) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
        //////////////////////////////////////////////////////////////////////////
        printf("\n vx mean max: %.1f at level: %d\n", vxmax, level);
    }
}
//////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridMeanLTwithDerivationsAndSqaredVelos(Parameter* para, int level, vector<string >& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("pressMed");
    nodedatanames.push_back("rhoMed");
    nodedatanames.push_back("vx1Med");
    nodedatanames.push_back("vx2Med");
    nodedatanames.push_back("vx3Med");
    nodedatanames.push_back("D(vx1)");
    nodedatanames.push_back("D(vx2)");
    nodedatanames.push_back("D(vx3)");
    nodedatanames.push_back("D(vx1)^2");
    nodedatanames.push_back("D(vx2)^2");
    nodedatanames.push_back("D(vx3)^2");
    nodedatanames.push_back("D(vx1*vx2)");
    nodedatanames.push_back("geo");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        nodedata[6].resize(sizeOfNodes);
        nodedata[7].resize(sizeOfNodes);
        nodedata[8].resize(sizeOfNodes);
        nodedata[9].resize(sizeOfNodes);
        nodedata[10].resize(sizeOfNodes);
        nodedata[11].resize(sizeOfNodes);
        nodedata[12].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");
        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1]  = para->getParH(level)->meanPressureOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[1][dn1]  = para->getParH(level)->meanDensityOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[2][dn1]  = para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio();
                nodedata[3][dn1]  = para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio();
                nodedata[4][dn1]  = para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio();
                nodedata[5][dn1]  = (((double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio()));
                nodedata[6][dn1]  = (((double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio()));
                nodedata[7][dn1]  = (((double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio()));
                nodedata[8][dn1]  = 
                    (((double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio())) * 
                    (((double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio()));
                nodedata[9][dn1]  = 
                    (((double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio())) * 
                    (((double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio()));
                nodedata[10][dn1] = 
                    (((double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio())) * 
                    (((double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio()));
                nodedata[11][dn1] =
                    (((double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio())) *
                    (((double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio()) - (para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio()));
                //nodedata[8][dn1]  = (((double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio()) * ((double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio()));
                //nodedata[9][dn1]  = (((double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio()) * ((double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio()));
                //nodedata[10][dn1] = (((double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio()) * ((double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio()));
                nodedata[12][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if (std::fabs(nodedata[2][dn1]) > std::fabs(vxmax)) vxmax = nodedata[2][dn1];
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid == true) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
        //////////////////////////////////////////////////////////////////////////
        printf("\n vx mean max: %.1f at level: %d\n", vxmax, level);
    }
}
//////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridEffMean(Parameter* para, int level, std::string& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    vector< string > nodedatanames;
    nodedatanames.push_back("pressMed");
    nodedatanames.push_back("rhoMed");
    nodedatanames.push_back("vx1Med");
    nodedatanames.push_back("vx2Med");
    nodedatanames.push_back("vx3Med");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    bool neighborsFluid;
    vector< vector<double>> nodedata(nodedatanames.size());

    nodes.resize(para->getParH(level)->numberOfNodes);
    nodedata[0].resize(para->getParH(level)->numberOfNodes);
    nodedata[1].resize(para->getParH(level)->numberOfNodes);
    nodedata[2].resize(para->getParH(level)->numberOfNodes);
    nodedata[3].resize(para->getParH(level)->numberOfNodes);
    nodedata[4].resize(para->getParH(level)->numberOfNodes);
    nodedata[5].resize(para->getParH(level)->numberOfNodes);

    for (size_t pos = 0; pos < para->getParH(level)->numberOfNodes; pos++)
    {
        if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
        {
            //////////////////////////////////////////////////////////////////////////
            double x1  = para->getParH(level)->coordinateX[pos];
            double x2  = para->getParH(level)->coordinateY[pos];
            double x3  = para->getParH(level)->coordinateZ[pos];
            double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
            double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
            double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
            //////////////////////////////////////////////////////////////////////////
            number1 = pos;
            neighborsFluid = true;
            //////////////////////////////////////////////////////////////////////////
            nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
            nodedata[0][number1] = para->getParH(level)->meanPressureOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[1][number1] = para->getParH(level)->meanDensityOut[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
            nodedata[2][number1] = para->getParH(level)->meanVelocityInXdirectionOut[pos] * para->getVelocityRatio();
            nodedata[3][number1] = para->getParH(level)->meanVelocityInYdirectionOut[pos] * para->getVelocityRatio();
            nodedata[4][number1] = para->getParH(level)->meanVelocityInZdirectionOut[pos] * para->getVelocityRatio();
            nodedata[5][number1] = para->getParH(level)->typeOfGridNode[pos];
            //////////////////////////////////////////////////////////////////////////
            number2 = para->getParH(level)->neighborX[number1];
            number3 = para->getParH(level)->neighborY[number2];
            number4 = para->getParH(level)->neighborY[number1];
            number5 = para->getParH(level)->neighborZ[number1];
            number6 = para->getParH(level)->neighborZ[number2];
            number7 = para->getParH(level)->neighborZ[number3];
            number8 = para->getParH(level)->neighborZ[number4];
            //////////////////////////////////////////////////////////////////////////
            if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            //if (level == 0 &&
            //    (number2 <= number1 ||
            //    number3 <= number1 ||
            //    number4 <= number1 ||
            //    number5 <= number1 ||
            //    number6 <= number1 ||
            //    number7 <= number1 ||
            //    number8 <= number1) )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
            //////////////////////////////////////////////////////////////////////////
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
}
//////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridEff2ndMoments(Parameter* para, int level, std::string& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    vector< string > nodedatanames;
    nodedatanames.push_back("kxyFromfcNEQ");
    nodedatanames.push_back("kyzFromfcNEQ");
    nodedatanames.push_back("kxzFromfcNEQ");
    nodedatanames.push_back("kxxMyyFromfcNEQ");
    nodedatanames.push_back("kxxMzzFromfcNEQ");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    bool neighborsFluid;
    vector< vector< double > > nodedata(nodedatanames.size());

    nodes.resize(para->getParH(level)->numberOfNodes);
    nodedata[0].resize(para->getParH(level)->numberOfNodes);
    nodedata[1].resize(para->getParH(level)->numberOfNodes);
    nodedata[2].resize(para->getParH(level)->numberOfNodes);
    nodedata[3].resize(para->getParH(level)->numberOfNodes);
    nodedata[4].resize(para->getParH(level)->numberOfNodes);
    nodedata[5].resize(para->getParH(level)->numberOfNodes);

    for (size_t pos = 0; pos < para->getParH(level)->numberOfNodes; pos++)
    {
        if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
        {
            //////////////////////////////////////////////////////////////////////////
            double x1  = para->getParH(level)->coordinateX[pos];
            double x2  = para->getParH(level)->coordinateY[pos];
            double x3  = para->getParH(level)->coordinateZ[pos];
            double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
            double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
            double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
            //////////////////////////////////////////////////////////////////////////
            number1 = pos;
            neighborsFluid = true;
            //////////////////////////////////////////////////////////////////////////
            nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
            nodedata[0][number1] = para->getParH(level)->kxyFromfcNEQ[pos];
            nodedata[1][number1] = para->getParH(level)->kyzFromfcNEQ[pos];
            nodedata[2][number1] = para->getParH(level)->kxzFromfcNEQ[pos];
            nodedata[3][number1] = para->getParH(level)->kxxMyyFromfcNEQ[pos];
            nodedata[4][number1] = para->getParH(level)->kxxMzzFromfcNEQ[pos];
            nodedata[5][number1] = para->getParH(level)->typeOfGridNode[pos];
            //////////////////////////////////////////////////////////////////////////
            number2 = para->getParH(level)->neighborX[number1];
            number3 = para->getParH(level)->neighborY[number2];
            number4 = para->getParH(level)->neighborY[number1];
            number5 = para->getParH(level)->neighborZ[number1];
            number6 = para->getParH(level)->neighborZ[number2];
            number7 = para->getParH(level)->neighborZ[number3];
            number8 = para->getParH(level)->neighborZ[number4];
            //////////////////////////////////////////////////////////////////////////
            if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
            //////////////////////////////////////////////////////////////////////////
            if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
            //////////////////////////////////////////////////////////////////////////
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname,nodes,cells,nodedatanames,nodedata);
}
//////////////////////////////////////////////////////////////////////////












//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridEff2ndMomentsLT(Parameter* para, int level, vector<string >& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("kxyFromfcNEQ");
    nodedatanames.push_back("kyzFromfcNEQ");
    nodedatanames.push_back("kxzFromfcNEQ");
    nodedatanames.push_back("kxxMyyFromfcNEQ");
    nodedatanames.push_back("kxxMzzFromfcNEQ");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part=0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        //printf("\n test in if II... \n");
        for (unsigned int pos=startpos;pos<endpos;pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
                double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
                double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][number1] = para->getParH(level)->kxyFromfcNEQ[pos];
                nodedata[1][number1] = para->getParH(level)->kyzFromfcNEQ[pos];
                nodedata[2][number1] = para->getParH(level)->kxzFromfcNEQ[pos];
                nodedata[3][number1] = para->getParH(level)->kxxMyyFromfcNEQ[pos];
                nodedata[4][number1] = para->getParH(level)->kxxMzzFromfcNEQ[pos];
                nodedata[5][number1] = para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part],nodes,cells,nodedatanames,nodedata);
    } 
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridEff3rdMomentsLT(Parameter* para, int level, vector<string >& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("CUMbbb");
    nodedatanames.push_back("CUMabc");
    nodedatanames.push_back("CUMbac");
    nodedatanames.push_back("CUMbca");
    nodedatanames.push_back("CUMcba");
    nodedatanames.push_back("CUMacb");
    nodedatanames.push_back("CUMcab");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part=0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        nodedata[6].resize(sizeOfNodes);
        nodedata[7].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");
        for (unsigned int pos=startpos;pos<endpos;pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
                double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
                double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][number1] = para->getParH(level)->CUMbbb[pos];
                nodedata[1][number1] = para->getParH(level)->CUMabc[pos];
                nodedata[2][number1] = para->getParH(level)->CUMbac[pos];
                nodedata[3][number1] = para->getParH(level)->CUMbca[pos];
                nodedata[4][number1] = para->getParH(level)->CUMcba[pos];
                nodedata[5][number1] = para->getParH(level)->CUMacb[pos];
                nodedata[6][number1] = para->getParH(level)->CUMcab[pos];
                nodedata[7][number1] = para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part],nodes,cells,nodedatanames,nodedata);
    } 
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeUnstrucuredGridEffHigherMomentsLT(Parameter* para, int level, vector<string >& fname) 
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleUInt8 > cells;
    //vector< UbTupleUInt8 > cells2;
    vector< string > nodedatanames;
    nodedatanames.push_back("CUMcbb");
    nodedatanames.push_back("CUMbcb");
    nodedatanames.push_back("CUMbbc");
    nodedatanames.push_back("CUMcca");
    nodedatanames.push_back("CUMcac");
    nodedatanames.push_back("CUMacc");
    nodedatanames.push_back("CUMbcc");
    nodedatanames.push_back("CUMcbc");
    nodedatanames.push_back("CUMccb");
    nodedatanames.push_back("CUMccc");
    nodedatanames.push_back("geo");
    unsigned int number1,number2,number3,number4,number5,number6,number7,number8;
    unsigned int dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8;
    bool neighborsFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    vector< vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part=0; part < fname.size(); part++)
    {
        vxmax = 0;
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        sizeOfNodes = FilePartCalculator::calculateNumberOfNodesInPart(para->getParH(level)->numberOfNodes);
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
        nodedata[6].resize(sizeOfNodes);
        nodedata[7].resize(sizeOfNodes);
        nodedata[8].resize(sizeOfNodes);
        nodedata[9].resize(sizeOfNodes);
        nodedata[10].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");
        for (unsigned int pos=startpos;pos<endpos;pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1  = para->getParH(level)->coordinateX[pos];
                double x2  = para->getParH(level)->coordinateY[pos];
                double x3  = para->getParH(level)->coordinateZ[pos];
                double x1P = para->getParH(level)->coordinateX[para->getParH(level)->neighborX[pos]];
                double x2P = para->getParH(level)->coordinateY[para->getParH(level)->neighborY[pos]];
                double x3P = para->getParH(level)->coordinateZ[para->getParH(level)->neighborZ[pos]];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[number1]=( makeUbTuple( (float)(x1 ),(float)(x2 ),(float)(x3 ) ) );
                nodedata[0][number1] = para->getParH(level)->CUMcbb[pos];
                nodedata[1][number1] = para->getParH(level)->CUMbcb[pos];
                nodedata[2][number1] = para->getParH(level)->CUMbbc[pos];
                nodedata[3][number1] = para->getParH(level)->CUMcca[pos];
                nodedata[4][number1] = para->getParH(level)->CUMcac[pos];
                nodedata[5][number1] = para->getParH(level)->CUMacc[pos];
                nodedata[6][number1] = para->getParH(level)->CUMbcc[pos];
                nodedata[7][number1] = para->getParH(level)->CUMcbc[pos];
                nodedata[8][number1] = para->getParH(level)->CUMccb[pos];
                nodedata[9][number1] = para->getParH(level)->CUMccc[pos];
                nodedata[10][number1] = para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID )  neighborsFluid=false;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid==true) cells.push_back( makeUbTuple(number1,number2,number3,number4,number5,number6,number7,number8) );        
                //////////////////////////////////////////////////////////////////////////
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part],nodes,cells,nodedatanames,nodedata);
    } 
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeQs(Parameter* para, int level, std::string& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleInt2 > qs;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    int node = 0;
    int wall = 1;
    int line = 0;
    double dx = 1.0 / pow(2, level);
    real* QQ;
    QforBoundaryConditions Q;
    double nodeX1, nodeX2, nodeX3, wallX1, wallX2, wallX3, q;
    //////////////////////////////////////////////////////////////////////////
    sizeOfNodes = para->getParH(level)->geometryBC.numberOfBCnodes;
    endpos = startpos + sizeOfNodes;
    //////////////////////////////////////////////////////////////////////////
    //qs.clear();
    //nodes.clear();
    unsigned int numberOfLines = para->getD3Qxx() * sizeOfNodes;
    unsigned int numberOfNodes = numberOfLines * 2;
    qs.resize(numberOfLines);
    nodes.resize(numberOfNodes);
    //////////////////////////////////////////////////////////////////////////
    vector< string > nodedatanames;
    nodedatanames.push_back("sizeQ");
    vector< vector< double > > nodedata(nodedatanames.size());
    nodedata[0].resize(numberOfNodes);
    //////////////////////////////////////////////////////////////////////////
    for (unsigned int pos = startpos; pos < endpos; pos++)
    {
        //////////////////////////////////////////////////////////////////////////
        nodeX1 = para->getParH(level)->coordinateX[para->getParH(level)->geometryBC.k[pos]];
        nodeX2 = para->getParH(level)->coordinateY[para->getParH(level)->geometryBC.k[pos]];
        nodeX3 = para->getParH(level)->coordinateZ[para->getParH(level)->geometryBC.k[pos]];
        wallX1 = 0.0;
        wallX2 = 0.0;
        wallX3 = 0.0;
        q      = 0.0;
        //////////////////////////////////////////////////////////////////////////
        for (size_t typeOfQ = vf::lbm::dir::STARTDIR; typeOfQ <= vf::lbm::dir::ENDDIR; typeOfQ++)
        {
            QQ = para->getParH(level)->geometryBC.q27[0];
            Q.q27[typeOfQ] = &QQ[typeOfQ*sizeOfNodes];
            q = (double)(Q.q27[typeOfQ][pos]);
            //////////////////////////////////////////////////////////////////////////
            switch (typeOfQ)
            {
                case E:   wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3;        break;
                case N:   wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case W:   wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3;        break;
                case S:   wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case NE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case NW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case SW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case SE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case T:   wallX1 = nodeX1;        wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TN:  wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case TW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TS:  wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case B:   wallX1 = nodeX1;        wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BN:  wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case BW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BS:  wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case TNE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case BSW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case BNE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case TSW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case TSE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case BNW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case BSE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case TNW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case REST:wallX1 = nodeX1;        wallX2 = nodeX2;         wallX3 = nodeX3;        break;
                default: throw UbException(UB_EXARGS, "unknown direction");
            }
            //////////////////////////////////////////////////////////////////////////
            nodes[node] = makeUbTuple((float)(nodeX1), (float)(nodeX2), (float)(nodeX3));
            nodes[wall] = makeUbTuple((float)(wallX1), (float)(wallX2), (float)(wallX3));
            qs[line]    = makeUbTuple(node, wall);
            //////////////////////////////////////////////////////////////////////////
            nodedata[0][node] = q;
            nodedata[0][wall] = q;
            //////////////////////////////////////////////////////////////////////////
            node = node + 2;
            wall = wall + 2;
            line++;
            //////////////////////////////////////////////////////////////////////////
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeLines(fname, nodes, qs);
    //WbWriterVtkXmlBinary::getInstance()->writeLinesWithNodeData(fname, nodes, qs, nodedatanames, nodedata);
    //WbWriterVtkXmlASCII::getInstance()->writeLinesWithNodeData(fname, nodes, qs, nodedatanames, nodedata);
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeQsInflow(Parameter* para, int level, std::string& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleInt2 > qs;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    int node = 0;
    int wall = 1;
    int line = 0;
    double dx = 1.0 / pow(2, level);
    real* QQ;
    QforBoundaryConditions Q;
    double nodeX1, nodeX2, nodeX3, wallX1, wallX2, wallX3, q;
    //////////////////////////////////////////////////////////////////////////
    sizeOfNodes = para->getParH(level)->velocityBC.numberOfBCnodes;
    endpos = startpos + sizeOfNodes;
    //////////////////////////////////////////////////////////////////////////
    //qs.clear();
    //nodes.clear();
    unsigned int numberOfLines = para->getD3Qxx() * sizeOfNodes;
    unsigned int numberOfNodes = numberOfLines * 2;
    qs.resize(numberOfLines);
    nodes.resize(numberOfNodes);
    //////////////////////////////////////////////////////////////////////////
    vector< string > nodedatanames;
    nodedatanames.push_back("sizeQ");
    vector< vector< double > > nodedata(nodedatanames.size());
    nodedata[0].resize(numberOfNodes);
    //////////////////////////////////////////////////////////////////////////
    for (unsigned int pos = startpos; pos < endpos; pos++)
    {
        //////////////////////////////////////////////////////////////////////////
        nodeX1 = para->getParH(level)->coordinateX[para->getParH(level)->velocityBC.k[pos]];
        nodeX2 = para->getParH(level)->coordinateY[para->getParH(level)->velocityBC.k[pos]];
        nodeX3 = para->getParH(level)->coordinateZ[para->getParH(level)->velocityBC.k[pos]];
        wallX1 = 0.0;
        wallX2 = 0.0;
        wallX3 = 0.0;
        q      = 0.0;
        //////////////////////////////////////////////////////////////////////////
        for (size_t typeOfQ = vf::lbm::dir::STARTDIR; typeOfQ <= vf::lbm::dir::ENDDIR; typeOfQ++)
        {
            QQ = para->getParH(level)->velocityBC.q27[0];
            Q.q27[typeOfQ] = &QQ[typeOfQ*sizeOfNodes];
            q = (double)(Q.q27[typeOfQ][pos]);
            if( q < 0.0 ) q = 0.0;
            //////////////////////////////////////////////////////////////////////////
            switch (typeOfQ)
            {
                case E:   wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3;        break;
                case N:   wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case W:   wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3;        break;
                case S:   wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case NE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case NW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case SW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case SE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case T:   wallX1 = nodeX1;        wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TN:  wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case TW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TS:  wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case B:   wallX1 = nodeX1;        wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BN:  wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case BW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BS:  wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case TNE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case BSW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case BNE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case TSW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case TSE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case BNW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case BSE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case TNW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case REST:wallX1 = nodeX1;        wallX2 = nodeX2;         wallX3 = nodeX3;        break;
                default: throw UbException(UB_EXARGS, "unknown direction");
            }
            //////////////////////////////////////////////////////////////////////////
            nodes[node] = makeUbTuple((float)(nodeX1), (float)(nodeX2), (float)(nodeX3));
            nodes[wall] = makeUbTuple((float)(wallX1), (float)(wallX2), (float)(wallX3));
            qs[line]    = makeUbTuple(node, wall);
            //////////////////////////////////////////////////////////////////////////
            nodedata[0][node] = q;
            nodedata[0][wall] = q;
            //////////////////////////////////////////////////////////////////////////
            node = node + 2;
            wall = wall + 2;
            line++;
            //////////////////////////////////////////////////////////////////////////
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeLines(fname, nodes, qs);
    //WbWriterVtkXmlBinary::getInstance()->writeLinesWithNodeData(fname, nodes, qs, nodedatanames, nodedata);
    //WbWriterVtkXmlASCII::getInstance()->writeLinesWithNodeData(fname, nodes, qs, nodedatanames, nodedata);
}
//////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////
void writeQsPressure(Parameter* para, int level, std::string& fname)
{
    vector< UbTupleFloat3 > nodes;
    vector< UbTupleInt2 > qs;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    int node = 0;
    int wall = 1;
    int line = 0;
    double dx = 1.0 / pow(2, level);
    real* QQ;
    QforBoundaryConditions Q;
    double nodeX1, nodeX2, nodeX3, wallX1, wallX2, wallX3, q;
    //////////////////////////////////////////////////////////////////////////
    sizeOfNodes = para->getParH(level)->pressureBC.numberOfBCnodes;
    endpos = startpos + sizeOfNodes;
    //////////////////////////////////////////////////////////////////////////
    //qs.clear();
    //nodes.clear();
    unsigned int numberOfLines = para->getD3Qxx() * sizeOfNodes;
    unsigned int numberOfNodes = numberOfLines * 2;
    qs.resize(numberOfLines);
    nodes.resize(numberOfNodes);
    //////////////////////////////////////////////////////////////////////////
    vector< string > nodedatanames;
    nodedatanames.push_back("sizeQ");
    vector< vector< double > > nodedata(nodedatanames.size());
    nodedata[0].resize(numberOfNodes);
    //////////////////////////////////////////////////////////////////////////
    for (unsigned int pos = startpos; pos < endpos; pos++)
    {
        //////////////////////////////////////////////////////////////////////////
        nodeX1 = para->getParH(level)->coordinateX[para->getParH(level)->pressureBC.k[pos]];
        nodeX2 = para->getParH(level)->coordinateY[para->getParH(level)->pressureBC.k[pos]];
        nodeX3 = para->getParH(level)->coordinateZ[para->getParH(level)->pressureBC.k[pos]];
        wallX1 = 0.0;
        wallX2 = 0.0;
        wallX3 = 0.0;
        q      = 0.0;
        //////////////////////////////////////////////////////////////////////////
        for (size_t typeOfQ = vf::lbm::dir::STARTDIR; typeOfQ <= vf::lbm::dir::ENDDIR; typeOfQ++)
        {
            QQ = para->getParH(level)->pressureBC.q27[0];
            Q.q27[typeOfQ] = &QQ[typeOfQ*sizeOfNodes];
            q = (double)(Q.q27[typeOfQ][pos]);
            if( q < 0.0 ) q = 0.0;
            //////////////////////////////////////////////////////////////////////////
            switch (typeOfQ)
            {
                case E:   wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3;        break;
                case N:   wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case W:   wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3;        break;
                case S:   wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case NE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case NW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3;        break;
                case SW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case SE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3;        break;
                case T:   wallX1 = nodeX1;        wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TN:  wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case TW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 + q*dx; break;
                case TS:  wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case B:   wallX1 = nodeX1;        wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BE:  wallX1 = nodeX1 + q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BN:  wallX1 = nodeX1;        wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case BW:  wallX1 = nodeX1 - q*dx; wallX2 = nodeX2;        wallX3 = nodeX3 - q*dx; break;
                case BS:  wallX1 = nodeX1;        wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case TNE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case BSW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case BNE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case TSW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case TSE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 + q*dx; break;
                case BNW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 - q*dx; break;
                case BSE: wallX1 = nodeX1 + q*dx; wallX2 = nodeX2 - q*dx; wallX3 = nodeX3 - q*dx; break;
                case TNW: wallX1 = nodeX1 - q*dx; wallX2 = nodeX2 + q*dx; wallX3 = nodeX3 + q*dx; break;
                case REST:wallX1 = nodeX1;        wallX2 = nodeX2;         wallX3 = nodeX3;        break;
                default: throw UbException(UB_EXARGS, "unknown direction");
            }
            //////////////////////////////////////////////////////////////////////////
            nodes[node] = makeUbTuple((float)(nodeX1), (float)(nodeX2), (float)(nodeX3));
            nodes[wall] = makeUbTuple((float)(wallX1), (float)(wallX2), (float)(wallX3));
            qs[line]    = makeUbTuple(node, wall);
            //////////////////////////////////////////////////////////////////////////
            nodedata[0][node] = q;
            nodedata[0][wall] = q;
            //////////////////////////////////////////////////////////////////////////
            node = node + 2;
            wall = wall + 2;
            line++;
            //////////////////////////////////////////////////////////////////////////
        }
    }
    WbWriterVtkXmlBinary::getInstance()->writeLines(fname, nodes, qs);
    //WbWriterVtkXmlBinary::getInstance()->writeLinesWithNodeData(fname, nodes, qs, nodedatanames, nodedata);
    //WbWriterVtkXmlASCII::getInstance()->writeLinesWithNodeData(fname, nodes, qs, nodedatanames, nodedata);
}
//////////////////////////////////////////////////////////////////////////
}

#endif
