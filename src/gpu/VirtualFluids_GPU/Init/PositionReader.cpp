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
		para->getParH(level)->QPropeller.numberOfBCnodes = in.readInteger();
		para->getParD(level)->QPropeller.numberOfBCnodes = para->getParH(level)->QPropeller.numberOfBCnodes;
		in.readLine();
		if (level == para->getFine())
		{
			for(int u=0; u<para->getParH(level)->QPropeller.numberOfBCnodes; u++)
			{
				test = in.readInteger();
				if (para->getParH(level)->geoSP[test] == GEO_FLUID)
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
			for(int u=0; u<para->getParH(level)->QPropeller.numberOfBCnodes; u++)
			{
				in.readInteger();
				in.readDouble();
				in.readDouble();
				in.readDouble();
				in.readLine();
			}
		}
		para->getParH(level)->QPropeller.numberOfBCnodes = count;
		para->getParD(level)->QPropeller.numberOfBCnodes = para->getParH(level)->QPropeller.numberOfBCnodes;
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
				if (para->getParH(level)->geoSP[test] == GEO_FLUID)
				{
					para->getParH(level)->QPropeller.k[count] = test; 
					para->getParH(level)->QPropeller.Vx[count] = (real)in.readDouble();
					para->getParH(level)->QPropeller.Vy[count] = (real)in.readDouble();
					para->getParH(level)->QPropeller.Vz[count] = (real)in.readDouble();
					para->getParH(level)->QPropeller.RhoBC[count] = 0.0f;									
					count++;
				}
				else
				{
					in.readDouble();
					in.readDouble();
					in.readDouble();
				}
				//para->getParH(level)->QPropeller.k[count] = test; 
				//para->getParH(level)->QPropeller.Vx[count] = (real)in.readDouble();
				//para->getParH(level)->QPropeller.Vy[count] = (real)in.readDouble();
				//para->getParH(level)->QPropeller.Vz[count] = (real)in.readDouble();
				//para->getParH(level)->QPropeller.Vx[count]	  = 0.07f;
				//para->getParH(level)->QPropeller.Vy[count]	  = 0.0f;
				//para->getParH(level)->QPropeller.Vz[count]	  = 0.0f;
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
	real* QQ                  = para->getParH(para->getFine())->QPropeller.q27[0]; 
	unsigned int sizeQ           = para->getParH(para->getFine())->QPropeller.numberOfBCnodes; 
	QforBoundaryConditions Q;
	Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
	Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
	Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
	Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
	Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
	Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
	Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
	Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
	Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
	Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
	Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
	Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
	Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
	Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
	Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
	Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
	Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
	Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
	Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
	Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
	Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
	Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
	Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
	Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
	Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
	Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
	Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
	//////////////////////////////////////////////////////////////////
	for(int u=0; u<para->getParH(para->getFine())->QPropeller.numberOfBCnodes; u++)
	{
		for (int dir = dirE; dir<=dirBSW; dir++)
		{
			if ((dir==dirE)  || 
				(dir==dirNE) || (dir==dirSE) || (dir==dirTE) || (dir==dirBE) ||
				(dir==dirTNE)|| (dir==dirBNE)|| (dir==dirTSE)|| (dir==dirBSE))
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
