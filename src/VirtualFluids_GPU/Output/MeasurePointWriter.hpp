#ifndef MEASURE_POINT_WRITER_H
#define MEASURE_POINT_WRITER_H

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

class MeasurePointWriter
{
public:
	MeasurePointWriter(){}
	~MeasurePointWriter(){}

	static void writeMeasurePoints(Parameter* para, int level, int index, int t)
	{
		ostringstream convert;   // stream used for the conversion
		convert << t;      // insert the textual representation of 'Number' in the characters in the stream
		std::string st = convert.str();
		UbFileOutputASCII out(para->getFName()+"_MeasurePoint_"+para->getParH(level)->MP[index].name+"_"+st+".dat");

		out.writeString("Level:");
		out.writeInteger(level);
		out.writeLine();
		out.writeString("Vx  Vy  Vz  Rho");
		out.writeLine();
		int numberNodes = (int)para->getParH(level)->MP[index].Rho.size();
		out.writeInteger(numberNodes);
		out.writeLine();
		for(int u=0; u<numberNodes; u++)
		{
			out.writeFloat((float)(para->getParH(level)->MP[index].Vx[u] * para->getVelocityRatio()));
			out.writeFloat((float)(para->getParH(level)->MP[index].Vy[u] * para->getVelocityRatio()));
			out.writeFloat((float)(para->getParH(level)->MP[index].Vz[u] * para->getVelocityRatio()));
			out.writeFloat((float)(para->getParH(level)->MP[index].Rho[u] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio()));
			out.writeLine();
		}
		out.writeLine();
	}

	static void writeTestAcousticXY(Parameter* para, int level, int t)
	{
		ostringstream convert;   // stream used for the conversion
		convert << t;      // insert the textual representation of 'Number' in the characters in the stream
		std::string st = convert.str();
		ostringstream convertLevel;   // stream used for the conversion
		convertLevel << level;      // insert the textual representation of 'Number' in the characters in the stream
		std::string sLevel = convertLevel.str();

		UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_XY.dat");

		int numberNodes = (int)para->getParH(level)->size_Mat_SP;
		
		doubflo deltaX = 1.0f / pow(2, level);
		doubflo halfDx = deltaX / 2.0f;
		doubflo middleOfTheGrid = (para->getMaxCoordZ()[0] + para->getMinCoordZ()[0]) / 2.f;
		//cout << "deltax: " << deltaX << ", halfDx: " << halfDx << ", middle of the grid: " << middleOfTheGrid << endl;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(level)->geoSP[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(level)->coordZ_SP[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(level)->coordZ_SP[u]) )
			{
				out.writeFloat((float)(para->getParH(level)->rho_SP[u]));
				out.writeFloat((float)(para->getParH(level)->press_SP[u]));
				out.writeFloat((float)(para->getParH(level)->coordX_SP[u]));
				out.writeFloat((float)(para->getParH(level)->coordY_SP[u]));
				out.writeLine();
			}
		}
	}

	static void writeTestAcousticYZ(Parameter* para, int level, int t)
	{
		ostringstream convert;   // stream used for the conversion
		convert << t;      // insert the textual representation of 'Number' in the characters in the stream
		std::string st = convert.str();
		ostringstream convertLevel;   // stream used for the conversion
		convertLevel << level;      // insert the textual representation of 'Number' in the characters in the stream
		std::string sLevel = convertLevel.str();

		UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_YZ.dat");

		int numberNodes = (int)para->getParH(level)->size_Mat_SP;

		doubflo deltaX = 1.0f / pow(2, level);
		doubflo halfDx = deltaX / 2.0f;
		doubflo middleOfTheGrid = (para->getMaxCoordX()[0] + para->getMinCoordX()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(level)->geoSP[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(level)->coordX_SP[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(level)->coordX_SP[u]))
			{
				out.writeFloat((float)(para->getParH(level)->rho_SP[u]));
				out.writeFloat((float)(para->getParH(level)->press_SP[u]));
				out.writeFloat((float)(para->getParH(level)->coordY_SP[u]));
				out.writeFloat((float)(para->getParH(level)->coordZ_SP[u]));
				out.writeLine();
			}
		}
	}

	static void writeTestAcousticXZ(Parameter* para, int level, int t)
	{
		ostringstream convert;   // stream used for the conversion
		convert << t;      // insert the textual representation of 'Number' in the characters in the stream
		std::string st = convert.str();
		ostringstream convertLevel;   // stream used for the conversion
		convertLevel << level;      // insert the textual representation of 'Number' in the characters in the stream
		std::string sLevel = convertLevel.str();

		UbFileOutputASCII out(para->getFName() + "_timestep_" + st + "_level_" + sLevel + "_XZ.dat");

		int numberNodes = (int)para->getParH(level)->size_Mat_SP;

		doubflo deltaX = 1.0f / pow(2, level);
		doubflo halfDx = deltaX / 2.0f;
		doubflo middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(level)->geoSP[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(level)->coordY_SP[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(level)->coordY_SP[u]))
			{
				out.writeFloat((float)(para->getParH(level)->rho_SP[u]));
				out.writeFloat((float)(para->getParH(level)->press_SP[u]));
				out.writeFloat((float)(para->getParH(level)->coordX_SP[u]));
				out.writeFloat((float)(para->getParH(level)->coordZ_SP[u]));
				out.writeLine();
			}
		}
	}


protected:

private:
};
#endif
