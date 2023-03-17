#ifndef POSITION_VECTOR_INTEGER_WRITER_H
#define POSITION_VECTOR_INTEGER_WRITER_H

//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
// #include <math.h>

//#include <cmath>

//#include "LBM/LB.h"
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
                out.writeInteger(para->getParH(level)->intCF.numberOfCells);
				out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->intCF.numberOfCells; u++)
				{
					out.writeInteger(para->getParH(level)->intCF.coarseCellIndices[u]);
				}
				out.writeLine();
			} //end levelloop
		}
		else if (Type == "_InterfaceCFF")
		{
			for (int level = 0; level < para->getMaxLevel(); level++)
			{
                out.writeInteger(para->getParH(level)->intCF.numberOfCells);
				out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->intCF.numberOfCells; u++)
				{
					out.writeInteger(para->getParH(level)->intCF.fineCellIndices[u]);
				}
				out.writeLine();
			} //end levelloop
		}
		else if (Type == "_InterfaceFCC")
		{
			for (int level = 0; level < para->getMaxLevel(); level++)
			{
                out.writeInteger(para->getParH(level)->intFC.numberOfCells);
				out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->intFC.numberOfCells; u++)
				{
					out.writeInteger(para->getParH(level)->intFC.coarseCellIndices[u]);
				}
				out.writeLine();
			} //end levelloop
		}
		else if (Type == "_InterfaceFCF")
		{
			for (int level = 0; level < para->getMaxLevel(); level++)
			{
                out.writeInteger(para->getParH(level)->intFC.numberOfCells);
				out.writeLine();
                for (unsigned int u = 0; u < para->getParH(level)->intFC.numberOfCells; u++)
				{
					out.writeInteger(para->getParH(level)->intFC.fineCellIndices[u]);
				}
				out.writeLine();
			} //end levelloop
		}
	}
protected:
private:
};
#endif
