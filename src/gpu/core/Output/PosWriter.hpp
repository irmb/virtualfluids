#ifndef POSITION_WRITER_H
#define POSITION_WRITER_H

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

class PositionWriter
{
public:
    PositionWriter(){}
    ~PositionWriter(){}

    static void writePosition(Parameter* para, std::string Type)
    {
        UbFileOutputASCII out(para->getFName()+Type+".dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        if (Type == "_geoSP")
        {
            for (int level = 0; level <= para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for(size_t index = 0; index < para->getParH(level)->numberOfNodes; index++)
                {
                    out.writeInteger(para->getParH(level)->typeOfGridNode[index]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_neighborX_SP")
        {
            for (int level = 0; level <= para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++)
                {
                    out.writeInteger(para->getParH(level)->neighborX[index]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_neighborY_SP")
        {
            for (int level = 0; level <= para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++)
                {
                    out.writeInteger(para->getParH(level)->neighborY[index]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_neighborZ_SP")
        {
            for (int level = 0; level <= para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->numberOfNodes);
                out.writeLine();
                for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++)
                {
                    out.writeInteger(para->getParH(level)->neighborZ[index]);
                }
                out.writeLine();
            } //end levelloop
        }
    }
protected:
private:
};
#endif
