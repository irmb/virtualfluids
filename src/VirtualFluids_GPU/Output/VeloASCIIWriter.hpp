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
		int numberNodes = (int)para->getParH(level)->size_Mat_SP;
		//write
		UbFileOutputASCII out(para->getFName() + "_VelocitiesASCII_" + std::to_string(level) + "_" + std::to_string(t) + ".dat");
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
			out.writeFloat((float)(para->getParH(level)->vx_SP[index] * para->getVelocityRatio()));
			out.writeFloat((float)(para->getParH(level)->vy_SP[index] * para->getVelocityRatio()));
			out.writeFloat((float)(para->getParH(level)->vz_SP[index] * para->getVelocityRatio()));
			out.writeLine();
		}
		out.writeLine();
	}

protected:

private:
};
#endif
