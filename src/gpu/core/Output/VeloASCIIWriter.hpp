#ifndef VELO_ASCII_WRITER_H
#define VELO_ASCII_WRITER_H

#include <numeric>
#include "basics/utilities/UbFileOutputASCII.h"
#include "Parameter/Parameter.h"

class VeloASCIIWriter
{
public:
    VeloASCIIWriter(){}
    ~VeloASCIIWriter(){}

    static void writeVelocitiesAsTXT(Parameter* para, int level, int t)
    {
        //calc
        int numberNodes = (int)para->getParH(level)->numberOfNodes;
        //write
        UbFileOutputASCII out(para->getFName() + "_VelocitiesASCII_" + std::to_string(level) + "_ID_" + StringUtil::toString<int>(para->getMyProcessID()) + "_t_" + std::to_string(t) + ".dat");
        //header
        out.writeString("Level:");
        out.writeInteger(level);
        out.writeLine();
        out.writeString("vX vY vZ");
        out.writeLine();
        //out.writeInteger(numberNodes);
        //out.writeLine();
        for (int index = 0; index < numberNodes; index++)
        {
            out.writeFloat((float)(para->getParH(level)->velocityX[index] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->velocityY[index] * para->getVelocityRatio()));
            out.writeFloat((float)(para->getParH(level)->velocityZ[index] * para->getVelocityRatio()));
            out.writeLine();
        }
        out.writeLine();
    }

protected:

private:
};
#endif
