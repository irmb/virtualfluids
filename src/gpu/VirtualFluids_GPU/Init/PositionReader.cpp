#include "PositionReader.h"

#include "Parameter/Parameter.h"

#include <basics/utilities/UbFileInputASCII.h>


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
	Q.q27[E   ] = &QQ[E   *sizeQ];
	Q.q27[W   ] = &QQ[W   *sizeQ];
	Q.q27[N   ] = &QQ[N   *sizeQ];
	Q.q27[S   ] = &QQ[S   *sizeQ];
	Q.q27[T   ] = &QQ[T   *sizeQ];
	Q.q27[B   ] = &QQ[B   *sizeQ];
	Q.q27[NE  ] = &QQ[NE  *sizeQ];
	Q.q27[SW  ] = &QQ[SW  *sizeQ];
	Q.q27[SE  ] = &QQ[SE  *sizeQ];
	Q.q27[NW  ] = &QQ[NW  *sizeQ];
	Q.q27[TE  ] = &QQ[TE  *sizeQ];
	Q.q27[BW  ] = &QQ[BW  *sizeQ];
	Q.q27[BE  ] = &QQ[BE  *sizeQ];
	Q.q27[TW  ] = &QQ[TW  *sizeQ];
	Q.q27[TN  ] = &QQ[TN  *sizeQ];
	Q.q27[BS  ] = &QQ[BS  *sizeQ];
	Q.q27[BN  ] = &QQ[BN  *sizeQ];
	Q.q27[TS  ] = &QQ[TS  *sizeQ];
	Q.q27[REST] = &QQ[REST*sizeQ];
	Q.q27[TNE ] = &QQ[TNE *sizeQ];
	Q.q27[TSW ] = &QQ[TSW *sizeQ];
	Q.q27[TSE ] = &QQ[TSE *sizeQ];
	Q.q27[TNW ] = &QQ[TNW *sizeQ];
	Q.q27[BNE ] = &QQ[BNE *sizeQ];
	Q.q27[BSW ] = &QQ[BSW *sizeQ];
	Q.q27[BSE ] = &QQ[BSE *sizeQ];
	Q.q27[BNW ] = &QQ[BNW *sizeQ];
	//////////////////////////////////////////////////////////////////
	for(uint u=0; u<para->getParH(para->getFine())->propellerBC.numberOfBCnodes; u++)
	{
		for (int dir = E; dir<=BSW; dir++)
		{
			if ((dir==E)  || 
				(dir==NE) || (dir==SE) || (dir==TE) || (dir==BE) ||
				(dir==TNE)|| (dir==BNE)|| (dir==TSE)|| (dir==BSE))
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
