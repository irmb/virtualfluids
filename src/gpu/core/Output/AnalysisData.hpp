#ifndef ANALYSIS_DATA_H
#define ANALYSIS_DATA_H

#include "basics/utilities/UbFileOutputASCII.h"
#include "Parameter/Parameter.h"
#include "StringUtilities/StringUtil.h"

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
				out.writeDouble((double)para->getParH(0)->velocityY[para->getParH(0)->k[m]]);
				//schr�g x
// 				out.writeDouble((double)para->getParH(0)->vz_SP[para->getParH(0)->k[m]]);
				//schr�g z
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
		int numberNodes = (int)para->getParH(0)->numberOfNodes;

		real deltaX = 1.0f;
		real halfDx = deltaX / 2.0f;
		real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(0)->typeOfGridNode[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(0)->coordinateY[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(0)->coordinateY[u]))
			{
				out.writeDouble((float)(para->getParH(0)->velocityX[u]));
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
				out.writeDouble((double)para->getParH(0)->velocityZ[para->getParH(0)->k[m]]);
			}
			out.writeLine();
		}
	}

	static void writeAnalysisDataXSP(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName() + "_AD_X_" + StringUtil::toString<int>(t) + ".dat");

		real level = 0; //uniform
		int numberNodes = (int)para->getParH(level)->numberOfNodes;

		real deltaX = 1.0f / pow(2, level);
		real halfDx = deltaX / 2.0f;
		real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u]))
			{
				out.writeDouble((double)(para->getParH(level)->velocityX[u]));
			}
		}
	}

	static void writeAnalysisDataYSP(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName() + "_AD_Y_" + StringUtil::toString<int>(t) + ".dat");

		real level = 0; //uniform
		int numberNodes = (int)para->getParH(level)->numberOfNodes;

		real deltaX = 1.0f / pow(2, level);
		real halfDx = deltaX / 2.0f;
		real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u]))
			{
				out.writeDouble((double)(para->getParH(level)->velocityY[u]));
			}
		}
	}

	static void writeAnalysisDataZSP(Parameter* para, unsigned int t)
	{
		UbFileOutputASCII out(para->getFName() + "_AD_Z_" + StringUtil::toString<int>(t) + ".dat");

		real level = 0; //uniform
		int numberNodes = (int)para->getParH(level)->numberOfNodes;

		real deltaX = 1.0f / pow(2, level);
		real halfDx = deltaX / 2.0f;
		real middleOfTheGrid = (para->getMaxCoordY()[0] + para->getMinCoordY()[0]) / 2.f;

		for (int u = 0; u < numberNodes; u++)
		{
			if ((para->getParH(level)->typeOfGridNode[u] == GEO_FLUID) &&
				((middleOfTheGrid - halfDx) <= para->getParH(level)->coordinateY[u]) &&
				((middleOfTheGrid + halfDx) >= para->getParH(level)->coordinateY[u]))
			{
				out.writeDouble((double)(para->getParH(level)->velocityZ[u]));
			}
		}
	}

protected:
private:
};
#endif
