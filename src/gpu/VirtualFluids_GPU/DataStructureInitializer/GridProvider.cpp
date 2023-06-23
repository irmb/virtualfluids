#include "GridProvider.h"

#include <Parameter/Parameter.h>
#include "GridReaderFiles/GridReader.h"
#include "GridReaderGenerator/GridGenerator.h"

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include <GPU/CudaMemoryManager.h>


std::shared_ptr<GridProvider> GridProvider::makeGridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager, vf::gpu::CommunicationRoutine& communicator)
{
    return std::shared_ptr<GridProvider>(new GridGenerator(builder, para, cudaMemoryManager, communicator));
}

std::shared_ptr<GridProvider> GridProvider::makeGridReader(FILEFORMAT format, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager)
{
    return std::shared_ptr<GridProvider>(new GridReader(format, para, cudaMemoryManager));
}

void GridProvider::setNumberOfNodes(uint numberOfNodes, int level) const
{
    para->getParH(level)->numberOfNodes          = (unsigned long long)numberOfNodes;
    para->getParD(level)->numberOfNodes          = (unsigned long long)numberOfNodes;
    para->getParH(level)->memSizeRealLBnodes     = sizeof(real) * para->getParH(level)->numberOfNodes;
    para->getParD(level)->memSizeRealLBnodes     = sizeof(real) * para->getParD(level)->numberOfNodes;
    para->getParH(level)->memSizeLonglongLBnodes = sizeof(unsigned long long) * para->getParH(level)->numberOfNodes;
    para->getParD(level)->memSizeLonglongLBnodes = sizeof(unsigned long long) * para->getParD(level)->numberOfNodes;
}

void GridProvider::setNumberOfTaggedFluidNodes(uint numberOfNodes, CollisionTemplate tag, int level) const
{
    para->getParH(level)->numberOfTaggedFluidNodes[tag] = numberOfNodes;
    para->getParD(level)->numberOfTaggedFluidNodes[tag] = numberOfNodes;
}

void GridProvider::setInitialNodeValues(uint numberOfNodes, int level) const
{
    for (uint pos = 1; pos <= numberOfNodes; pos++)
    {
        const real coordX = para->getParH(level)->coordinateX[pos];
        const real coordY = para->getParH(level)->coordinateY[pos];
        const real coordZ = para->getParH(level)->coordinateZ[pos];

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

        para->getParH(level)->rho[pos] = rho; 
        para->getParH(level)->velocityX[pos]  = vx; 
        para->getParH(level)->velocityY[pos]  = vy;
        para->getParH(level)->velocityZ[pos]  = vz; 

        //////////////////////////////////////////////////////////////////////////

        if (para->getCalcMedian()) {
            para->getParH(level)->vx_SP_Med[pos] = 0.0f;
            para->getParH(level)->vy_SP_Med[pos] = 0.0f;
            para->getParH(level)->vz_SP_Med[pos] = 0.0f;
            para->getParH(level)->rho_SP_Med[pos] = 0.0f;
            para->getParH(level)->press_SP_Med[pos] = 0.0f;
        }
        if (para->getUseWale()) {
            para->getParH(level)->turbViscosity[pos] = 0.0f;
            //Debug
            para->getParH(level)->gSij[pos] = 0.0f;
            para->getParH(level)->gSDij[pos] = 0.0f;
            para->getParH(level)->gDxvx[pos] = 0.0f;
            para->getParH(level)->gDyvx[pos] = 0.0f;
            para->getParH(level)->gDzvx[pos] = 0.0f;
            para->getParH(level)->gDxvy[pos] = 0.0f;
            para->getParH(level)->gDyvy[pos] = 0.0f;
            para->getParH(level)->gDzvy[pos] = 0.0f;
            para->getParH(level)->gDxvz[pos] = 0.0f;
            para->getParH(level)->gDyvz[pos] = 0.0f;
            para->getParH(level)->gDzvz[pos] = 0.0f;
        }

        if (para->getIsBodyForce()) {
            para->getParH(level)->forceX_SP[pos] = 0.0f;
            para->getParH(level)->forceY_SP[pos] = 0.0f;
            para->getParH(level)->forceZ_SP[pos] = 0.0f;
        }
    }


}


void GridProvider::setPressSizePerLevel(int level, int sizePerLevel) const
{
    para->getParH(level)->pressureBC.numberOfBCnodes = sizePerLevel;
    para->getParD(level)->pressureBC.numberOfBCnodes = sizePerLevel;
    para->getParH(level)->numberOfPressureBCnodesRead = sizePerLevel * para->getD3Qxx();
    para->getParD(level)->numberOfPressureBCnodesRead = sizePerLevel * para->getD3Qxx();
}


void GridProvider::setVelocitySizePerLevel(int level, int sizePerLevel) const
{
    para->getParH(level)->velocityBC.numberOfBCnodes = sizePerLevel;
    para->getParD(level)->velocityBC.numberOfBCnodes = sizePerLevel;
    para->getParH(level)->numberOfVeloBCnodesRead = sizePerLevel * para->getD3Qxx();
    para->getParD(level)->numberOfVeloBCnodesRead = sizePerLevel * para->getD3Qxx();
}

void GridProvider::setOutflowSizePerLevel(int level, int sizePerLevel) const
{
    para->getParH(level)->outflowBC.numberOfBCnodes = sizePerLevel;
    para->getParD(level)->outflowBC.numberOfBCnodes = sizePerLevel;
    para->getParH(level)->numberOfOutflowBCnodesRead = sizePerLevel * para->getD3Qxx();
    para->getParD(level)->numberOfOutflowBCnodesRead = sizePerLevel * para->getD3Qxx();
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
