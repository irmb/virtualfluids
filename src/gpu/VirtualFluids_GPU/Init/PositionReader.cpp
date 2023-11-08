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
	Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
	Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
	Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
	Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
	Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
	Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
	Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
	Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
	Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
	Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
	Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
	Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
	Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
	Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
	Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
	Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
	Q.q27[d000] = &QQ[d000*sizeQ];
	Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
	Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
	Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
	Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
	Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
	Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
	Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
	Q.q27[dMPM ] = &QQ[dMPM *sizeQ];
	//////////////////////////////////////////////////////////////////
	for(uint u=0; u<para->getParH(para->getFine())->propellerBC.numberOfBCnodes; u++)
	{
		for (size_t dir = dP00; dir<=dMMM; dir++)
		{
			if ((dir==dP00)  || 
				(dir==dPP0) || (dir==dPM0) || (dir==dP0P) || (dir==dP0M) ||
				(dir==dPPP)|| (dir==dPPM)|| (dir==dPMP)|| (dir==dPMM))
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
