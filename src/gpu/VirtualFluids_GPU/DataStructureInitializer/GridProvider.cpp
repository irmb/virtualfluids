#include "GridProvider.h"

#include <Parameter/Parameter.h>
#include "GridReaderFiles/GridReader.h"
#include "GridReaderGenerator/GridGenerator.h"

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include <GPU/CudaMemoryManager.h>


std::shared_ptr<GridProvider> GridProvider::makeGridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager)
{
    return std::shared_ptr<GridProvider>(new GridGenerator(builder, para, cudaManager));
}

std::shared_ptr<GridProvider> GridProvider::makeGridReader(FILEFORMAT format, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager)
{
    return std::shared_ptr<GridProvider>(new GridReader(format, para, cudaManager));
}

void GridProvider::setNumberOfNodes(const int numberOfNodes, const int level) const
{
    para->getParH(level)->size_Mat_SP = numberOfNodes;
    para->getParD(level)->size_Mat_SP = numberOfNodes;
    para->getParH(level)->mem_size_real_SP = sizeof(real) * para->getParH(level)->size_Mat_SP;
    para->getParH(level)->mem_size_int_SP = sizeof(uint) * para->getParH(level)->size_Mat_SP;
    para->getParD(level)->mem_size_real_SP = sizeof(real) * para->getParD(level)->size_Mat_SP;
    para->getParD(level)->mem_size_int_SP = sizeof(uint) * para->getParD(level)->size_Mat_SP;
}

void GridProvider::setInitalNodeValues(const int numberOfNodes, const int level) const
{
    for (int j = 1; j <= numberOfNodes; j++)
    {
        const real coordX = para->getParH(level)->coordX_SP[j];
        const real coordY = para->getParH(level)->coordY_SP[j];
        const real coordZ = para->getParH(level)->coordZ_SP[j];

        real rho, vx, vy, vz;

        // call functor object with initial condition
        if( para->getInitialCondition() )
        {
            para->getInitialCondition()(coordX,coordY,coordZ,rho,vx,vy,vz);
        }
        else
        {
            rho = real(0.0);
            vx  = real(0.0);
            vy  = real(0.0);
            vz  = real(0.0);
        }

        para->getParH(level)->rho_SP[j] = rho; 
        para->getParH(level)->vx_SP[j]  = vx; 
        para->getParH(level)->vy_SP[j]  = vy;
        para->getParH(level)->vz_SP[j]  = vz; 

        //////////////////////////////////////////////////////////////////////////

        if (para->getCalcMedian()) {
            para->getParH(level)->vx_SP_Med[j] = 0.0f;
            para->getParH(level)->vy_SP_Med[j] = 0.0f;
            para->getParH(level)->vz_SP_Med[j] = 0.0f;
            para->getParH(level)->rho_SP_Med[j] = 0.0f;
            para->getParH(level)->press_SP_Med[j] = 0.0f;
        }
        if (para->getUseWale()) {
            para->getParH(level)->turbViscosity[j] = 0.0f;
            //Debug
            para->getParH(level)->gSij[j] = 0.0f;
            para->getParH(level)->gSDij[j] = 0.0f;
            para->getParH(level)->gDxvx[j] = 0.0f;
            para->getParH(level)->gDyvx[j] = 0.0f;
            para->getParH(level)->gDzvx[j] = 0.0f;
            para->getParH(level)->gDxvy[j] = 0.0f;
            para->getParH(level)->gDyvy[j] = 0.0f;
            para->getParH(level)->gDzvy[j] = 0.0f;
            para->getParH(level)->gDxvz[j] = 0.0f;
            para->getParH(level)->gDyvz[j] = 0.0f;
            para->getParH(level)->gDzvz[j] = 0.0f;
        }
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

    for (int level = para->getCoarse(); level <= para->getFine(); level++)
    {
        cudaMemoryManager->cudaAllocLevelForcing(level);
        cudaMemoryManager->cudaCopyLevelForcingToDevice(level);
    }
}

void GridProvider::allocAndCopyQuadricLimiters()
{
    cudaMemoryManager->cudaAllocQuadricLimiters();
    cudaMemoryManager->cudaCopyQuadricLimitersToDevice();
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
