#include "PositionReader.h"

#include "Parameter/Parameter.h"

#include <basics/utilities/UbFileInputASCII.h>

using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePropellerCylinderForAlloc(Parameter* para)
{
	UbFileInputASCII in(para->getpropellerCylinder());
	int test = 0, count = 0;
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		para->getParH(level)->propellerBC.numberOfBCnodes = in.readInteger();
		para->getParD(level)->propellerBC.numberOfBCnodes = para->getParH(level)->propellerBC.numberOfBCnodes;
		in.readLine();
		if (level == para->getFine())
		{
			for(uint u=0; u<para->getParH(level)->propellerBC.numberOfBCnodes; u++)
			{
				test = in.readInteger();
				if (para->getParH(level)->typeOfGridNode[test] == GEO_FLUID)
				{
					count++;
				}
				////////////////////////////////////////////////////////////////////////
				//for(unsigned int ix3=0; ix3<para->getParH(level)->nz; ix3++)
				//{
				//	for(unsigned int ix2=0; ix2<para->getParH(level)->ny; ix2++)
				//	{
				//		for(unsigned int ix1=0; ix1<para->getParH(level)->nx; ix1++)
				//		{
				//			unsigned int m = para->getParH(level)->nx*(para->getParH(level)->ny*ix3 + ix2) + ix1;
				//			if (para->getParH(level)->k[m] == test)
				//			{
				//				if(para->getParH(level)->geo[m] == 1)
				//				{
				//					count++;									
				//				}
				//			}
				//		}
				//	}
				//}
				//count++;
				////////////////////////////////////////////////////////////////////////
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		else
		{
			for(uint u=0; u<para->getParH(level)->propellerBC.numberOfBCnodes; u++)
			{
				in.readInteger();
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		para->getParH(level)->propellerBC.numberOfBCnodes = count;
		para->getParD(level)->propellerBC.numberOfBCnodes = para->getParH(level)->propellerBC.numberOfBCnodes;
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readFilePropellerCylinder(Parameter* para)
{
	UbFileInputASCII in(para->getpropellerCylinder());
	int test = 0, count = 0;
	int maxlevel = in.readInteger();
	in.readLine();

	for (int level = 0; level < maxlevel; level++)
	{
		int allnodes = in.readInteger();
		in.readLine();
		if (level == para->getFine())
		{
			for(int u=0; u<allnodes; u++)
			{
				test = in.readInteger();
				////////////////////////////////////////////////////////////////////////
				if (para->getParH(level)->typeOfGridNode[test] == GEO_FLUID)
				{
					para->getParH(level)->propellerBC.k[count] = test; 
					para->getParH(level)->propellerBC.Vx[count] = (real)in.readDouble();
					para->getParH(level)->propellerBC.Vy[count] = (real)in.readDouble();
					para->getParH(level)->propellerBC.Vz[count] = (real)in.readDouble();
					para->getParH(level)->propellerBC.RhoBC[count] = 0.0f;									
					count++;
				}
				else
				{
					in.readDouble();
					in.readDouble();
					in.readDouble();
				}
				//para->getParH(level)->propellerBC.k[count] = test; 
				//para->getParH(level)->propellerBC.Vx[count] = (real)in.readDouble();
				//para->getParH(level)->propellerBC.Vy[count] = (real)in.readDouble();
				//para->getParH(level)->propellerBC.Vz[count] = (real)in.readDouble();
				//para->getParH(level)->propellerBC.Vx[count]	  = 0.07f;
				//para->getParH(level)->propellerBC.Vy[count]	  = 0.0f;
				//para->getParH(level)->propellerBC.Vz[count]	  = 0.0f;
				in.readLine();
			}
		} 
		else
		{
			for(int u=0; u<allnodes; u++)
			{
				in.readInteger(); 
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		printf("allnodes = %d, count = %d\n", allnodes, count);
	}
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::definePropellerQs(Parameter* para)
{
	//////////////////////////////////////////////////////////////////
	//preprocessing
	real* QQ                  = para->getParH(para->getFine())->propellerBC.q27[0]; 
	unsigned int sizeQ           = para->getParH(para->getFine())->propellerBC.numberOfBCnodes; 
	QforBoundaryConditions Q;
	Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
	Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
	Q.q27[DIR_0P0   ] = &QQ[DIR_0P0   *sizeQ];
	Q.q27[DIR_0M0   ] = &QQ[DIR_0M0   *sizeQ];
	Q.q27[DIR_00P   ] = &QQ[DIR_00P   *sizeQ];
	Q.q27[DIR_00M   ] = &QQ[DIR_00M   *sizeQ];
	Q.q27[DIR_PP0  ] = &QQ[DIR_PP0  *sizeQ];
	Q.q27[DIR_MM0  ] = &QQ[DIR_MM0  *sizeQ];
	Q.q27[DIR_PM0  ] = &QQ[DIR_PM0  *sizeQ];
	Q.q27[DIR_MP0  ] = &QQ[DIR_MP0  *sizeQ];
	Q.q27[DIR_P0P  ] = &QQ[DIR_P0P  *sizeQ];
	Q.q27[DIR_M0M  ] = &QQ[DIR_M0M  *sizeQ];
	Q.q27[DIR_P0M  ] = &QQ[DIR_P0M  *sizeQ];
	Q.q27[DIR_M0P  ] = &QQ[DIR_M0P  *sizeQ];
	Q.q27[DIR_0PP  ] = &QQ[DIR_0PP  *sizeQ];
	Q.q27[DIR_0MM  ] = &QQ[DIR_0MM  *sizeQ];
	Q.q27[DIR_0PM  ] = &QQ[DIR_0PM  *sizeQ];
	Q.q27[DIR_0MP  ] = &QQ[DIR_0MP  *sizeQ];
	Q.q27[d000] = &QQ[d000*sizeQ];
	Q.q27[DIR_PPP ] = &QQ[DIR_PPP *sizeQ];
	Q.q27[DIR_MMP ] = &QQ[DIR_MMP *sizeQ];
	Q.q27[DIR_PMP ] = &QQ[DIR_PMP *sizeQ];
	Q.q27[DIR_MPP ] = &QQ[DIR_MPP *sizeQ];
	Q.q27[DIR_PPM ] = &QQ[DIR_PPM *sizeQ];
	Q.q27[DIR_MMM ] = &QQ[DIR_MMM *sizeQ];
	Q.q27[DIR_PMM ] = &QQ[DIR_PMM *sizeQ];
	Q.q27[DIR_MPM ] = &QQ[DIR_MPM *sizeQ];
	//////////////////////////////////////////////////////////////////
	for(uint u=0; u<para->getParH(para->getFine())->propellerBC.numberOfBCnodes; u++)
	{
		for (size_t dir = dP00; dir<=DIR_MMM; dir++)
		{
			if ((dir==dP00)  || 
				(dir==DIR_PP0) || (dir==DIR_PM0) || (dir==DIR_P0P) || (dir==DIR_P0M) ||
				(dir==DIR_PPP)|| (dir==DIR_PPM)|| (dir==DIR_PMP)|| (dir==DIR_PMM))
			{
				Q.q27[dir][u] = 1.0f;
			} 
			else
			{
				Q.q27[dir][u] = -1.0f;
			}
		}
	}
	//////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
void PositionReader::readMeasurePoints( Parameter* para ) 
{
	UbFileInputASCII in(para->getmeasurePoints());
	int numberOfAllNodes = in.readInteger();
	in.readLine();
	int tempLevel;
	MeasurePoints tempMP;
	//printf("done, init the values...\n");
	for (int u = 0; u < numberOfAllNodes; u++)
	{
		tempMP.name = in.readString(); 		
		//printf("done, read the name...\n");
		tempMP.k = in.readInteger();
		//printf("done, read k...\n");
		tempLevel = in.readInteger();
		//printf("done, read level...\n");
		in.readLine();
		//printf("done, read the values...\n");
		para->getParH(tempLevel)->MP.push_back(tempMP);
		//printf("done, put it into a vector...\n");
	}
}
