#ifndef POSITION_WRITER_H
#define POSITION_WRITER_H

//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
// #include <math.h>

//#include <cmath>

//#include "LBM/LB.h"
//#include "LBM/D3Q27.h"
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
				for(unsigned int u=0; u<para->getParH(level)->numberOfNodes; u++)
				{
					out.writeInteger(para->getParH(level)->typeOfGridNode[u]);
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
				for(unsigned int u=0; u<para->getParH(level)->numberOfNodes; u++)
				{
					out.writeInteger(para->getParH(level)->neighborX[u]);
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
				for(unsigned int u=0; u<para->getParH(level)->numberOfNodes; u++)
				{
					out.writeInteger(para->getParH(level)->neighborY[u]);
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
				for(unsigned int u=0; u<para->getParH(level)->numberOfNodes; u++)
				{
					out.writeInteger(para->getParH(level)->neighborZ[u]);
				}
				out.writeLine();
			} //end levelloop
		}
	}
protected:
private:
};
#endif
