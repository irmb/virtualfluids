#ifndef ANALYSIS_DATA_H
#define ANALYSIS_DATA_H

#include "basics/utilities/UbFileOutputASCII.h"
#include "Parameter/Parameter.h"
#include "Utilities/StringUtil.hpp"

class AnalysisData
{
public:
	AnalysisData(){}
	~AnalysisData(){}

	static void writeAnalysisData(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName()+"_AD_"+StringUtil::toString<int>(t)+".dat");

		unsigned int j = para->getParH(0)->gridNY / 2 + STARTOFFY;

		for (unsigned int k=STARTOFFZ; k<para->getParH(0)->gridNZ + STARTOFFZ; k++)
		{
			for (unsigned int i=STARTOFFX; i<para->getParH(0)->gridNX + STARTOFFX; i++)
			{
				int m = para->getParH(0)->nx*(para->getParH(0)->ny*k + j) + i;
				//gerade
				out.writeDouble((double)para->getParH(0)->vy_SP[para->getParH(0)->k[m]]);
				//schräg x
// 				out.writeDouble((double)para->getParH(0)->vz_SP[para->getParH(0)->k[m]]);
				//schräg z
				//out.writeDouble((double)para->getParH(0)->vx_SP[para->getParH(0)->k[m]]);
			}
			out.writeLine();
		}
	}

	static void writeAnalysisDataX(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName()+"_AD_X_"+StringUtil::toString<int>(t)+".dat");

		//unsigned int j = para->getParH(0)->gridNY / 2 + STARTOFFY;

		//for (unsigned int k=STARTOFFZ; k<para->getParH(0)->gridNZ + STARTOFFZ; k++)
		//{
		//	for (unsigned int i=STARTOFFX; i<para->getParH(0)->gridNX + STARTOFFX; i++)
		//	{
		//		int m = para->getParH(0)->nx*(para->getParH(0)->ny*k + j) + i;
		//		//Taylor Green Vortex - X
		//		out.writeDouble((double)para->getParH(0)->vx_SP[para->getParH(0)->k[m]]);
		//	}
		//	out.writeLine();
		//}
		int numberNodes = (int)para->getParH(0)->size_Mat_SP;

		doubflo deltaX = 1.0f;
		doubflo halfDx = deltaX / 2.0f;
		doubflo middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(0)->geoSP[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(0)->coordY_SP[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(0)->coordY_SP[u]))
			{
				out.writeDouble((float)(para->getParH(0)->vx_SP[u]));
				out.writeLine();
			}
		}
	}

	static void writeAnalysisDataZ(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName()+"_AD_Z_"+StringUtil::toString<int>(t)+".dat");

		unsigned int j = para->getParH(0)->gridNY / 2 + STARTOFFY;

		for (unsigned int k=STARTOFFZ; k<para->getParH(0)->gridNZ + STARTOFFZ; k++)
		{
			for (unsigned int i=STARTOFFX; i<para->getParH(0)->gridNX + STARTOFFX; i++)
			{
				int m = para->getParH(0)->nx*(para->getParH(0)->ny*k + j) + i;
				//Taylor Green Vortex - Z
				out.writeDouble((double)para->getParH(0)->vz_SP[para->getParH(0)->k[m]]);
			}
			out.writeLine();
		}
	}

	static void writeAnalysisDataXSP(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName() + "_AD_X_" + StringUtil::toString<int>(t) + ".dat");

		doubflo level = 0; //uniform
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
				out.writeDouble((double)(para->getParH(level)->vx_SP[u]));
			}
		}
	}

	static void writeAnalysisDataYSP(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName() + "_AD_Y_" + StringUtil::toString<int>(t) + ".dat");

		doubflo level = 0; //uniform
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
				out.writeDouble((double)(para->getParH(level)->vy_SP[u]));
			}
		}
	}

	static void writeAnalysisDataZSP(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName() + "_AD_Z_" + StringUtil::toString<int>(t) + ".dat");

		doubflo level = 0; //uniform
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
				out.writeDouble((double)(para->getParH(level)->vz_SP[u]));
			}
		}
	}

protected:
private:
};
#endif
