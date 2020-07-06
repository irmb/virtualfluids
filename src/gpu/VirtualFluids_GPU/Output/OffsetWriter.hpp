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
				out.writeInteger(para->getParH(level)->K_CF);
				out.writeLine();
				for(unsigned int u=0; u<para->getParH(level)->K_CF; u++)
				{
					out.writeDouble(para->getParH(level)->offCF.xOffCF[u]);
					out.writeDouble(para->getParH(level)->offCF.yOffCF[u]);
					out.writeDouble(para->getParH(level)->offCF.zOffCF[u]);
				}
				out.writeLine();
			} //end levelloop
		}
		else if (Type == "_OffsetFC")
		{
			for (int level = 0; level < para->getMaxLevel(); level++)
			{
				out.writeInteger(para->getParH(level)->K_FC);
				out.writeLine();
				for(unsigned int u=0; u<para->getParH(level)->K_FC; u++)
				{
					out.writeDouble(para->getParH(level)->offFC.xOffFC[u]);
					out.writeDouble(para->getParH(level)->offFC.yOffFC[u]);
					out.writeDouble(para->getParH(level)->offFC.zOffFC[u]);
				}
				out.writeLine();
			} //end levelloop
		}
	}
protected:
private:
};
#endif
