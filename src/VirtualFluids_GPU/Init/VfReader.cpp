#include "Init/VfReader.h"

////////////////////////////////////////////////////////////////////////////////
void readVFkFull(Parameter* para, const std::string geometryFile)
{
	kFullReader::readFileForAlloc(geometryFile, para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaAllocFull(lev);
		//////////////////////////////////////////////////////////////////////////
		for(unsigned int ix3=0; ix3<para->getParH(lev)->nz; ix3++)
		{
			for(unsigned int ix2=0; ix2<para->getParH(lev)->ny; ix2++)
			{
				for(unsigned int ix1=0; ix1<para->getParH(lev)->nx; ix1++)
				{
					unsigned int m = para->getParH(lev)->nx*(para->getParH(lev)->ny*ix3 + ix2) + ix1;
					para->getParH(lev)->k[m]   =  0;
					para->getParH(lev)->geo[m] =  16;
				}
			}
		}
	}

	kFullReader::readFile(geometryFile, para);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readVFgeoFull(Parameter* para, const std::string geometryFile)
{
	kFullReader::readGeoFull(geometryFile, para);

	//////////////////////////////////////////////////////////////////////////
	//for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	//{
	//	for(unsigned int ix3=0; ix3<para->getParH(lev)->nz; ix3++)
	//	{
	//		for(unsigned int ix2=0; ix2<para->getParH(lev)->ny; ix2++)
	//		{
	//			for(unsigned int ix1=0; ix1<para->getParH(lev)->nx; ix1++)
	//			{
	//				unsigned int m = para->getParH(lev)->nx*(para->getParH(lev)->ny*ix3 + ix2) + ix1;
	//				if (para->getParH(lev)->geo[m] == 0 || para->getParH(lev)->geo[m] == 15)
	//				{
	//					para->getParH(lev)->geo[m] = 16;
	//				}
	//			}
	//		}
	//	}
	//}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readVecSP(Parameter* para)
{
	PositionReader::readFileForAlloc(para->getgeoVec(), para);

	//alloc
	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaAllocSP(lev);
	}

	//geoSP
	PositionReader::readFile(para->getgeoVec(), "geoVec", para);
	//neighborX
	PositionReader::readFile(para->getneighborX(), "neighborX", para);
	//neighborY
	PositionReader::readFile(para->getneighborY(), "neighborY", para);
	//neighborZ
	PositionReader::readFile(para->getneighborZ(), "neighborZ", para);

	//Copy Host -> Device
	for(int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		for(unsigned int u=0; u<para->getParH(lev)->size_Mat_SP; u++)
		{
			para->getParH(lev)->rho_SP[u]   = 0.01f;//+ lev/100.f;
			para->getParH(lev)->vx_SP[u]    = 0.0f;//+ lev/100.f;   
			para->getParH(lev)->vy_SP[u]    = 0.0f;//+ lev/100.f;   
			para->getParH(lev)->vz_SP[u]    = 0.0f;//+ lev/100.f;   
			para->getParH(lev)->press_SP[u] = 0.0f;//+ lev/100.f;
			//Median
			para->getParH(lev)->rho_SP_Med[u]   = 0.0f;
			para->getParH(lev)->vx_SP_Med[u]    = 0.0f;
			para->getParH(lev)->vy_SP_Med[u]    = 0.0f;
			para->getParH(lev)->vz_SP_Med[u]    = 0.0f;
			para->getParH(lev)->press_SP_Med[u] = 0.0f;
		}
		para->cudaCopySP(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readInterfaceCF(Parameter* para)
{
	PositionReader::readFileInterfaceForAlloc(para->getscaleCFC(), "CF", para);

	//alloc
	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaAllocInterfaceCF(lev);
	}

	//Scale Coarse to Fine - Coarse
	PositionReader::readFileInterface(para->getscaleCFC(), "CFC", para);
	//Scale Coarse to Fine - Fine
	PositionReader::readFileInterface(para->getscaleCFF(), "CFF", para);

	//Copy Host -> Device
	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaCopyInterfaceCF(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readInterfaceFC(Parameter* para)
{
	PositionReader::readFileInterfaceForAlloc(para->getscaleFCC(), "FC", para);

	//alloc
	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaAllocInterfaceFC(lev);
	}

	//Scale Fine to Coarse - Coarse
	PositionReader::readFileInterface(para->getscaleFCC(), "FCC", para);
	//Scale Fine to Coarse - Fine
	PositionReader::readFileInterface(para->getscaleFCF(), "FCF", para);

	//Copy Host -> Device
	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaCopyInterfaceFC(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readInterfaceOffCF(Parameter* para, const std::string geometryFile)
{
	PositionReader::readFileInterfaceOffsetForAlloc(geometryFile, "CF", para);

	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaAllocInterfaceOffCF(lev);
	}

	PositionReader::readFileInterfaceOffset(geometryFile, "CF", para);

	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaCopyInterfaceOffCF(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readInterfaceOffFC(Parameter* para, const std::string geometryFile)
{
	PositionReader::readFileInterfaceOffsetForAlloc(geometryFile, "FC", para);

	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaAllocInterfaceOffFC(lev);
	}

	PositionReader::readFileInterfaceOffset(geometryFile, "FC", para);

	for (int lev = 0; lev < para->getMaxLevel(); lev++)
	{
		para->cudaCopyInterfaceOffFC(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readNoSlipBc(Parameter* para)
{
	PositionReader::readFileNoSlipBcForAlloc(para->getnoSlipBcPos(), para);
	PositionReader::readFileNoSlipBcQreadForAlloc(para->getnoSlipBcQs(), para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaAllocWallBC(lev);
	}

	PositionReader::readFileNoSlipBcPos(para->getnoSlipBcPos(), para);
	PositionReader::readFileNoSlipBcValue(para->getnoSlipBcValue(), para);
	PositionReader::readFileNoSlipBcQs(para->getnoSlipBcQs(), para);

	PositionReader::findQs(para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaCopyWallBC(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readSlipBc(Parameter* para)
{
	PositionReader::readFileSlipBcForAlloc(para->getslipBcPos(), para);
	PositionReader::readFileSlipBcQreadForAlloc(para->getslipBcQs(), para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaAllocSlipBC(lev);
	}

	PositionReader::readFileSlipBcPos(para->getslipBcPos(), para);
	PositionReader::readFileSlipBcValue(para->getslipBcValue(), para);
	PositionReader::readFileSlipBcQs(para->getslipBcQs(), para);

	PositionReader::findSlipQs(para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaCopySlipBC(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readPressBc(Parameter* para)
{
	PositionReader::readFilePressBcForAlloc(para->getpressBcPos(), para);
	PositionReader::readFilePressBcQreadForAlloc(para->getpressBcQs(), para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaAllocPress(lev);
	}
	//only Coarse
	//para->cudaAllocPress(para->getCoarse());

	PositionReader::readFilePressBcPos(para->getpressBcPos(), para);
	PositionReader::readFilePressBcValue(para->getpressBcValue(), para);
	PositionReader::readFilePressBcQs(para->getpressBcQs(), para);

	PositionReader::findPressQs(para);

	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		para->cudaCopyPress(lev);
	}
	//only Coarse
	//para->cudaCopyPress(para->getCoarse());
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readPropellerCylinder(Parameter* para)
{
	PositionReader::readFilePropellerCylinderForAlloc(para);

	para->cudaAllocVeloPropeller(para->getFine());

	PositionReader::readFilePropellerCylinder(para);
	//PositionReader::definePropellerQs(para);

	para->cudaCopyVeloPropeller(para->getFine());
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void readMeasurePoints(Parameter* para)
{
	//read measure points from file
	PositionReader::readMeasurePoints(para);
	//printf("done, reading the file...\n");
	//level loop
	for (int lev = 0; lev <= para->getMaxLevel(); lev++)
	{
		//set Memory Size and malloc of the indices and macroscopic values per level
		para->getParH(lev)->numberOfValuesMP = (unsigned int)para->getParH(lev)->MP.size()*(unsigned int)para->getclockCycleForMP()/((unsigned int)para->getTimestepForMP());
		para->getParD(lev)->numberOfValuesMP = para->getParH(lev)->numberOfValuesMP;

		para->getParH(lev)->numberOfPointskMP = (int)para->getParH(lev)->MP.size();
		para->getParD(lev)->numberOfPointskMP = para->getParH(lev)->numberOfPointskMP;

		para->getParH(lev)->memSizeIntkMP = sizeof(unsigned int)*(int)para->getParH(lev)->MP.size();
		para->getParD(lev)->memSizeIntkMP = para->getParH(lev)->memSizeIntkMP;

		para->getParH(lev)->memSizeDoubflokMP = sizeof(doubflo)*para->getParH(lev)->numberOfValuesMP;
		para->getParD(lev)->memSizeDoubflokMP = para->getParH(lev)->memSizeDoubflokMP;		
		
		printf("Level: %d, numberOfValuesMP: %d, memSizeIntkMP: %d, memSizeDoubflokMP: %d\n",lev,para->getParH(lev)->numberOfValuesMP,para->getParH(lev)->memSizeIntkMP, para->getParD(lev)->memSizeDoubflokMP);

		para->cudaAllocMeasurePointsIndex(lev);

		//loop over all measure points per level 
		for(int index = 0; index < (int)para->getParH(lev)->MP.size(); index++)
		{
			//set indices
			para->getParH(lev)->kMP[index] = para->getParH(lev)->MP[index].k;
		}
		//loop over all measure points per level times MPClockCycle
		for(int index = 0; index < (int)para->getParH(lev)->numberOfValuesMP; index++)
		{
			//init values
			para->getParH(lev)->VxMP[index]  = (doubflo)0.0;
			para->getParH(lev)->VyMP[index]  = (doubflo)0.0;
			para->getParH(lev)->VzMP[index]  = (doubflo)0.0;
			para->getParH(lev)->RhoMP[index] = (doubflo)0.0;
		}

		//copy indices-arrays
		para->cudaCopyMeasurePointsIndex(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////








