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
                out.writeInteger(para->getParH(level)->intCF.kCF);
				out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->intCF.kCF; u++)
				{
					out.writeDouble(para->getParH(level)->neighborCF.x[u]);
					out.writeDouble(para->getParH(level)->neighborCF.y[u]);
					out.writeDouble(para->getParH(level)->neighborCF.z[u]);
				}
				out.writeLine();
			} //end levelloop
		}
		else if (Type == "_OffsetFC")
		{
			for (int level = 0; level < para->getMaxLevel(); level++)
			{
                out.writeInteger(para->getParH(level)->intFC.kFC);
				out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->intFC.kFC; u++)
				{
					out.writeDouble(para->getParH(level)->neighborFC.x[u]);
					out.writeDouble(para->getParH(level)->neighborFC.y[u]);
					out.writeDouble(para->getParH(level)->neighborFC.z[u]);
				}
				out.writeLine();
			} //end levelloop
		}
	}
protected:
private:
};
#endif
