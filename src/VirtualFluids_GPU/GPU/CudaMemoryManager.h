#ifndef CudamemoryManager_H
#define CudamemoryManager_H

#include <vector>
#include <string>
#include <memory>

#include "LBM/LB.h"
#include "LBM/D3Q27.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

class Parameter;

class CudaMemoryManager
{
public:
	static std::shared_ptr<CudaMemoryManager> make(std::shared_ptr<Parameter> parameter);
	
	void cudaAllocCoord(int lev);
	void cudaCopyCoord(int lev);
	void cudaFreeCoord(int lev);

	void cudaCopyDataToHost(int lev);

	void cudaAllocSP(int lev);
	void cudaCopySP(int lev);
	void cudaFreeSP(int lev);

	void cudaAllocVeloBC(int lev);
	void cudaCopyVeloBC(int lev);
	void cudaFreeVeloBC(int lev);
	void cudaAllocOutflowBC(int lev);
	void cudaCopyOutflowBC(int lev);
	void cudaFreeOutflowBC(int lev);
	void cudaAllocWallBC(int lev);
	void cudaCopyWallBC(int lev);
	void cudaFreeWallBC(int lev);

	void cudaAllocGeomBC(int lev);
	void cudaCopyGeomBC(int lev);
	void cudaFreeGeomBC(int lev);

	void cudaAllocPress(int lev);
	void cudaCopyPress(int lev);
	void cudaFreePress(int lev);

	void cudaAllocForcing();
	void cudaCopyForcingToDevice();
	void cudaCopyForcingToHost();
	void cudaFreeForcing();

	//////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	void cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborX(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborY(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborYIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborY(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborZ(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborZIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor);
	//////////////////////////////////////////////////////////////////////////

//
//private:
//	int coarse, fine, maxlevel;
//	int factor_gridNZ;
//	int D3Qxx;
//	InitCondition ic;
//	double memsizeGPU;
//	unsigned int limitOfNodesForVTK;

private:
    CudaMemoryManager(std::shared_ptr<Parameter> parameter);
    CudaMemoryManager(const CudaMemoryManager&);

    std::shared_ptr<Parameter> parameter;

};

#endif

