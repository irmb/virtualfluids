#ifndef POSITION_VECTOR_INTEGER_WRITER_H
#define POSITION_VECTOR_INTEGER_WRITER_H

//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
// #include <math.h>

//#include <cmath>

//#include "Calculation/Calculation.h"
//#include "lbm/constants/D3Q27.h"
#include "basics/utilities/UbFileOutputASCII.h"
#include "Parameter/Parameter.h"

//using namespace std;

//namespace kFullWriter
//{
//
//
//}
class PositionVectorIntegerWriter
{
public:
    PositionVectorIntegerWriter(){}
    ~PositionVectorIntegerWriter(){}

    static void writePositionInterface(Parameter* para, std::string Type)
    {
        UbFileOutputASCII out(para->getFName()+Type+".dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        if (Type == "_InterfaceCFC")
        {
            for (int level = 0; level < para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->coarseToFine.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++)
                {
                    out.writeInteger(para->getParH(level)->coarseToFine.coarseCellIndices[u]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_InterfaceCFF")
        {
            for (int level = 0; level < para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->coarseToFine.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++)
                {
                    out.writeInteger(para->getParH(level)->coarseToFine.fineCellIndices[u]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_InterfaceFCC")
        {
            for (int level = 0; level < para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->fineToCoarse.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++)
                {
                    out.writeInteger(para->getParH(level)->fineToCoarse.coarseCellIndices[u]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_InterfaceFCF")
        {
            for (int level = 0; level < para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->fineToCoarse.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++)
                {
                    out.writeInteger(para->getParH(level)->fineToCoarse.fineCellIndices[u]);
                }
                out.writeLine();
            } //end levelloop
        }
    }
protected:
private:
};
#endif
