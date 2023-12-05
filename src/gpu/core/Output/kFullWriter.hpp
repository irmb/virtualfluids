#ifndef K_FULL_WRITER_H
#define K_FULL_WRITER_H

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
class kFullWriter
{
public:
    kFullWriter(){}
    ~kFullWriter(){}

    static void writekFull(Parameter* para)
    {
        UbFileOutputASCII out(para->getFName()+"_kfull.dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        for (int level = 0; level <= para->getMaxLevel(); level++)
        {
            int nodeNumberX1 = para->getParH(level)->gridNX;
            int nodeNumberX2 = para->getParH(level)->gridNY;
            int nodeNumberX3 = para->getParH(level)->gridNZ;
            out.writeInteger(nodeNumberX1);
            out.writeInteger(nodeNumberX2);
            out.writeInteger(nodeNumberX3);
            out.writeLine();
            for(int ix3=0; ix3<nodeNumberX3; ix3++)
            {
                for(int ix2=0; ix2<nodeNumberX2; ix2++)
                {
                    for(int ix1=0; ix1<nodeNumberX1; ix1++)
                    {
                        unsigned int m = nodeNumberX1*(nodeNumberX2*ix3 + ix2) + ix1;
                        int number = para->getParH(level)->k[m];

                        //if(number<0)  out.writeInteger(0); 
                        //else          out.writeInteger(number+1);
                        out.writeInteger(number);
                    }
                }
            }
            out.writeLine();
        } //end levelloop

    }
protected:
private:
};
#endif
