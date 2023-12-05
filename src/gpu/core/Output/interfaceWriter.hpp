#ifndef INTERFACE_WRITER_H
#define INTERFACE_WRITER_H

//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
// #include <math.h>

//#include <cmath>

//#include "Calculation/Calculation.h"
//#include "lbm/constants/D3Q27.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbSystem.h>
#include "Parameter/Parameter.h"

//das stimmt alles noch nicht....
class interfaceWriter
{
public:
    interfaceWriter(){}
    ~interfaceWriter(){}

    static void writeVtkInterfaceCells(Parameter* para, std::string Type)
    {
        std::vector< UbTupleFloat3 > nodesVec;
        std::vector< UbTupleInt8 > cellsVec;
        int nodeNumberVec = 0;
        for (int level = 0; level < para->getMaxLevel(); level++)
        {
            if ((Type == "_InterfaceCFC") || (Type == "_InterfaceCFF"))
            {
                nodeNumberVec += para->getParH(level)->coarseToFine.numberOfCells;
            }
            else if (Type == "_InterfaceFCF")
            {
                nodeNumberVec += para->getParH(level)->fineToCoarse.numberOfCells;
            }
        }
        nodesVec.resize(nodeNumberVec*8);
        int nodeCount = 0;
        for (int level = 0; level < para->getMaxLevel(); level++)
        {
            int nx1lev = para->getParH(level)->gridNX;//((grid->getNX1()<<(level+levelDiff))*blocknx1)+1;
            int nx2lev = para->getParH(level)->gridNY;//((grid->getNX2()<<(level+levelDiff))*blocknx2)+1;
            int nx3lev = para->getParH(level)->gridNZ;//((grid->getNX3()<<(level+levelDiff))*blocknx3)+1;

            //int nodeDifferenz = nodeCounter;
            double nodeDeltaLevel = para->getParH(level)->dx; //nodeDelta/(1<<(level+levelDiff));
            double halfNodeDeltaLevel = 0.5*nodeDeltaLevel;
            //int count = 0;

            //std::vector<unsigned int>& posVec = posIndexVec[level];
            if (Type == "_InterfaceCFC")
            {
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++)
                {
                    int pos = para->getParH(level)->coarseToFine.coarseCellIndices[u];
                    int ix1 = pos % nx1lev;
                    int wertDurchNx1 = pos / nx1lev;
                    int ix2 = wertDurchNx1 % nx2lev;
                    int ix3 = wertDurchNx1 / nx2lev;

                    //cout<<"danach:"<<ix1<<" "<<ix2<<" "<<ix3<<" -> "<<pos<<endl;
                    double x1 = ix1*nodeDeltaLevel+halfNodeDeltaLevel;
                    double x2 = ix2*nodeDeltaLevel+halfNodeDeltaLevel;
                    double x3 = ix3*nodeDeltaLevel+halfNodeDeltaLevel;
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2+nodeDeltaLevel),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+nodeDeltaLevel),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2+nodeDeltaLevel),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+nodeDeltaLevel),(float)(x3+nodeDeltaLevel) ) );

                    cellsVec.push_back( makeUbTuple(nodeCount-8,nodeCount-7,nodeCount-6,nodeCount-5,nodeCount-4,nodeCount-3,nodeCount-2,nodeCount-1) );

                }                
            }
            else if (Type == "_InterfaceCFF")
            {
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++)
                {
                    int pos = para->getParH(level)->coarseToFine.fineCellIndices[u];
                    int ix1 = pos % nx1lev;
                    int wertDurchNx1 = pos / nx1lev;
                    int ix2 = wertDurchNx1 % nx2lev;
                    int ix3 = wertDurchNx1 / nx2lev;

                    //cout<<"danach:"<<ix1<<" "<<ix2<<" "<<ix3<<" -> "<<pos<<endl;
                    double x1 = ix1*nodeDeltaLevel+halfNodeDeltaLevel;
                    double x2 = ix2*nodeDeltaLevel+halfNodeDeltaLevel;
                    double x3 = ix3*nodeDeltaLevel+halfNodeDeltaLevel;
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2+nodeDeltaLevel),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+nodeDeltaLevel),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2+nodeDeltaLevel),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+nodeDeltaLevel),(float)(x3+nodeDeltaLevel) ) );

                    cellsVec.push_back( makeUbTuple(nodeCount-8,nodeCount-7,nodeCount-6,nodeCount-5,nodeCount-4,nodeCount-3,nodeCount-2,nodeCount-1) );

                }                
            }
            else if (Type == "_InterfaceFCF")
            {
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++)
                {
                    int pos = para->getParH(level)->fineToCoarse.fineCellIndices[u];
                    int ix1 = pos % nx1lev;
                    int wertDurchNx1 = pos / nx1lev;
                    int ix2 = wertDurchNx1 % nx2lev;
                    int ix3 = wertDurchNx1 / nx2lev;

                    //cout<<"danach:"<<ix1<<" "<<ix2<<" "<<ix3<<" -> "<<pos<<endl;
                    double x1 = ix1*nodeDeltaLevel+halfNodeDeltaLevel;
                    double x2 = ix2*nodeDeltaLevel+halfNodeDeltaLevel;
                    double x3 = ix3*nodeDeltaLevel+halfNodeDeltaLevel;
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2+nodeDeltaLevel),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+nodeDeltaLevel),(float)(x3) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+nodeDeltaLevel),(float)(x2+nodeDeltaLevel),(float)(x3+nodeDeltaLevel) ) );
                    nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+nodeDeltaLevel),(float)(x3+nodeDeltaLevel) ) );

                    cellsVec.push_back( makeUbTuple(nodeCount-8,nodeCount-7,nodeCount-6,nodeCount-5,nodeCount-4,nodeCount-3,nodeCount-2,nodeCount-1) );
                }
            }
        }
        // WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(filename,nodes,cells,nodedatanames,nodedata);
        std::string filenameVec = para->getFName()+Type+".dat";
        WbWriterVtkXmlBinary::getInstance()->writeOcts(filenameVec,nodesVec,cellsVec);
    }

    static void writeVtkInterfaceNodes(Parameter* para, std::string Type)
    {
        std::vector< UbTupleFloat3 > nodesVec;
        std::vector< UbTupleInt8 > cellsVec;
        int nodeNumberVec = 0;
        for (int level = 0; level < para->getMaxLevel(); level++)
        {
            nodeNumberVec += para->getParH(level)->fineToCoarse.numberOfCells;
        }
        nodesVec.resize(nodeNumberVec*8);
        int nodeCount = 0;
        for (int level = 0; level < para->getMaxLevel(); level++)
        {
            int nx1lev =  para->getParH(level)->gridNX;//((grid->getNX1()<<(level+levelDiff))*blocknx1)+1;
            int nx2lev =  para->getParH(level)->gridNY;//((grid->getNX2()<<(level+levelDiff))*blocknx2)+1;
            int nx3lev =  para->getParH(level)->gridNZ;//((grid->getNX3()<<(level+levelDiff))*blocknx3)+1;

            //int nodeDifferenz = nodeCounter;
            double nodeDeltaLevel = para->getParH(level)->dx; //nodeDelta/(1<<(level+levelDiff));
            double halfNodeDeltaLevel = 0.5*nodeDeltaLevel;
            double viertelNodeDelta = 0.25*nodeDeltaLevel;
            double achtelNodeDelta = 0.125*nodeDeltaLevel;
            //int count = 0;
            //std::vector<unsigned int>& posVec = posIndexVec[level];
            for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++)
            {
                int pos = para->getParH(level)->fineToCoarse.coarseCellIndices[u];//posVec[u];
                int ix1 = pos % nx1lev;
                int wertDurchNx1 = pos / nx1lev;
                int ix2 = wertDurchNx1 % nx2lev;
                int ix3 = wertDurchNx1 / nx2lev;

                //cout<<"danach:"<<ix1<<" "<<ix2<<" "<<ix3<<" -> "<<pos<<endl;
                double x1 = ix1*nodeDeltaLevel+halfNodeDeltaLevel-achtelNodeDelta;
                double x2 = ix2*nodeDeltaLevel+halfNodeDeltaLevel-achtelNodeDelta;
                double x3 = ix3*nodeDeltaLevel+halfNodeDeltaLevel-achtelNodeDelta;
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+viertelNodeDelta),(float)(x2),(float)(x3) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+viertelNodeDelta),(float)(x2+viertelNodeDelta),(float)(x3) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+viertelNodeDelta),(float)(x3) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2),(float)(x3+viertelNodeDelta) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+viertelNodeDelta),(float)(x2),(float)(x3+viertelNodeDelta) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1+viertelNodeDelta),(float)(x2+viertelNodeDelta),(float)(x3+viertelNodeDelta) ) );
                nodesVec[nodeCount++]=( makeUbTuple( (float)(x1),(float)(x2+viertelNodeDelta),(float)(x3+viertelNodeDelta) ) );

                cellsVec.push_back( makeUbTuple(nodeCount-8,nodeCount-7,nodeCount-6,nodeCount-5,nodeCount-4,nodeCount-3,nodeCount-2,nodeCount-1) );

            }
        }
        // WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(filename,nodes,cells,nodedatanames,nodedata);
        std::string filenameVec = para->getFName()+Type+".dat";
        WbWriterVtkXmlBinary::getInstance()->writeOcts(filenameVec,nodesVec,cellsVec);
    }

protected:
private:
};
#endif
