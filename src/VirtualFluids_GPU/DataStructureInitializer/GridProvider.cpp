#include "GridProvider.h"

#include <Parameter/Parameter.h>
#include "GridReaderFiles/GridReader.h"
#include "GridReaderGenerator/GridGenerator.h"

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include <GPU/CudaMemoryManager.h>


std::shared_ptr<GridProvider> GridProvider::makeGridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<GridProvider>(new GridGenerator(builder, para));
}

std::shared_ptr<GridProvider> GridProvider::makeGridReader(bool readBinaryFiles, std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<GridProvider>(new GridReader(readBinaryFiles, para));
}

void GridProvider::setNumberOfNodes(const int numberOfNodes, const int level) const
{
	para->getParH(level)->size_Mat_SP = numberOfNodes;
	para->getParD(level)->size_Mat_SP = numberOfNodes;
	para->getParH(level)->mem_size_doubflo_SP = sizeof(doubflo) * para->getParH(level)->size_Mat_SP;
	para->getParH(level)->mem_size_int_SP = sizeof(unsigned int) * para->getParH(level)->size_Mat_SP;
	para->getParD(level)->mem_size_doubflo_SP = sizeof(doubflo) * para->getParD(level)->size_Mat_SP;
	para->getParD(level)->mem_size_int_SP = sizeof(unsigned int) * para->getParD(level)->size_Mat_SP;
}

void GridProvider::setInitalNodeValues(const int numberOfNodes, const int level) const
{
	for (int j = 0; j <= numberOfNodes; j++)
	{
		para->getParH(level)->vx_SP[j] = para->getVelocity();//0.0f;//0.035f;
		para->getParH(level)->vy_SP[j] = 0.0f;//para->getVelocity();//0.0f;
		para->getParH(level)->vz_SP[j] = 0.0f;
		para->getParH(level)->rho_SP[j] = 0.0f;
		para->getParH(level)->press_SP[j] = 0.0f;
	}
}


void GridProvider::setPressSizePerLevel(int level, int sizePerLevel) const
{
	para->getParH(level)->QPress.kQ = sizePerLevel;
	para->getParD(level)->QPress.kQ = sizePerLevel;
	para->getParH(level)->kPressQread = sizePerLevel * para->getD3Qxx();
	para->getParD(level)->kPressQread = sizePerLevel * para->getD3Qxx();
}


void GridProvider::setVelocitySizePerLevel(int level, int sizePerLevel) const
{
	para->getParH(level)->Qinflow.kQ = sizePerLevel;
	para->getParD(level)->Qinflow.kQ = sizePerLevel;
	para->getParH(level)->kInflowQ = sizePerLevel;
	para->getParD(level)->kInflowQ = sizePerLevel;
	para->getParH(level)->kInflowQread = sizePerLevel * para->getD3Qxx();
	para->getParD(level)->kInflowQread = sizePerLevel * para->getD3Qxx();
}

void GridProvider::setOutflowSizePerLevel(int level, int sizePerLevel) const
{
	para->getParH(level)->Qoutflow.kQ = sizePerLevel;
	para->getParD(level)->Qoutflow.kQ = sizePerLevel;
	para->getParH(level)->kOutflowQread = sizePerLevel * para->getD3Qxx();
	para->getParD(level)->kOutflowQread = sizePerLevel * para->getD3Qxx();
}

void GridProvider::allocAndCopyForcing()
{
    cudaMemoryManager->cudaAllocForcing();
    cudaMemoryManager->cudaCopyForcingToDevice();
}

void GridProvider::freeMemoryOnHost()
{
    for (int level = para->getCoarse(); level <= para->getFine(); level++)
    {
        cudaMemoryManager->cudaFreeCoord(level);
        cudaMemoryManager->cudaFreeSP(level);
    }
}

void GridProvider::cudaCopyDataToHost(int level)
{
    cudaMemoryManager->cudaCopyDataToHost(level);
}
