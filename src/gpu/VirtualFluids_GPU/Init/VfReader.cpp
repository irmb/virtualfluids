#include "Init/VfReader.h"

#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "Init/PositionReader.h"
#include "GPU/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////
void readPropellerCylinder(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
	PositionReader::readFilePropellerCylinderForAlloc(para);

	cudaMemoryManager->cudaAllocVeloPropeller(para->getFine());

	PositionReader::readFilePropellerCylinder(para);
	//PositionReader::definePropellerQs(para);

	cudaMemoryManager->cudaCopyVeloPropeller(para->getFine());
}

////////////////////////////////////////////////////////////////////////////////
void readMeasurePoints(Parameter* para, CudaMemoryManager* cudaMemoryManager)
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

		para->getParH(lev)->memSizerealkMP = sizeof(real)*para->getParH(lev)->numberOfValuesMP;
		para->getParD(lev)->memSizerealkMP = para->getParH(lev)->memSizerealkMP;		
		
		printf("Level: %d, numberOfValuesMP: %d, memSizeIntkMP: %d, memSizerealkMP: %d\n",lev,para->getParH(lev)->numberOfValuesMP,para->getParH(lev)->memSizeIntkMP, para->getParD(lev)->memSizerealkMP);

		cudaMemoryManager->cudaAllocMeasurePointsIndex(lev);

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
			para->getParH(lev)->VxMP[index]  = (real)0.0;
			para->getParH(lev)->VyMP[index]  = (real)0.0;
			para->getParH(lev)->VzMP[index]  = (real)0.0;
			para->getParH(lev)->RhoMP[index] = (real)0.0;
		}

		//copy indices-arrays
		cudaMemoryManager->cudaCopyMeasurePointsIndex(lev);
	}
}
////////////////////////////////////////////////////////////////////////////////








