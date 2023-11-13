#ifndef OFFSET_WRITER_H
#define OFFSET_WRITER_H

#include "basics/utilities/UbFileOutputASCII.h"
#include "Parameter/Parameter.h"

class OffsetWriter
{
public:
    OffsetWriter(){}
    ~OffsetWriter(){}

    static void writeOffset(Parameter* para, std::string Type)
    {
        UbFileOutputASCII out(para->getFName()+Type+".dat");

        out.writeInteger(para->getMaxLevel());
        out.writeLine();

        if (Type == "_OffsetCF")
        {
            for (int level = 0; level < para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->coarseToFine.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->coarseToFine.numberOfCells; u++)
                {
                    out.writeDouble(para->getParH(level)->neighborCoarseToFine.x[u]);
                    out.writeDouble(para->getParH(level)->neighborCoarseToFine.y[u]);
                    out.writeDouble(para->getParH(level)->neighborCoarseToFine.z[u]);
                }
                out.writeLine();
            } //end levelloop
        }
        else if (Type == "_OffsetFC")
        {
            for (int level = 0; level < para->getMaxLevel(); level++)
            {
                out.writeInteger(para->getParH(level)->fineToCoarse.numberOfCells);
                out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->fineToCoarse.numberOfCells; u++)
                {
                    out.writeDouble(para->getParH(level)->neighborFineToCoarse.x[u]);
                    out.writeDouble(para->getParH(level)->neighborFineToCoarse.y[u]);
                    out.writeDouble(para->getParH(level)->neighborFineToCoarse.z[u]);
                }
                out.writeLine();
            } //end levelloop
        }
    }
protected:
private:
};
#endif
