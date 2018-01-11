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
    para->getParH(level)->mem_size_real_SP = sizeof(real) * para->getParH(level)->size_Mat_SP;
    para->getParH(level)->mem_size_int_SP = sizeof(unsigned int) * para->getParH(level)->size_Mat_SP;
    para->getParD(level)->mem_size_real_SP = sizeof(real) * para->getParD(level)->size_Mat_SP;
    para->getParD(level)->mem_size_int_SP = sizeof(unsigned int) * para->getParD(level)->size_Mat_SP;
}

void GridProvider::setInitalNodeValues(const int numberOfNodes, const int level) const
{
    //Taylor Green Vortex uniform
    real PI = 3.141592653589793238462643383279f;
    real gridX = para->getParH(level)->gridNX - 1;
    real gridY = para->getParH(level)->gridNY - 1;
    real gridZ = para->getParH(level)->gridNZ - 1;
    //like MG
    //real uAdvect = real (1. / 250.); //32 nodes -> 250; 40 nodes -> 200; 64 nodes -> 500; 128 nodes -> 1000; 256 nodes -> 2000; 512 nodes -> 4000
    real uAdvect = real(0.0016); //32 nodes -> 0.032; 64 nodes -> 0.016; 128 nodes -> 0.008; 256 nodes -> 0.004; 512 nodes -> 0.002

    for (int j = 0; j <= numberOfNodes; j++)
    {
        para->getParH(level)->rho_SP[j] =
            (real)((para->getVelocity()*para->getVelocity())*3.0 / 4.0*(cos(para->getParH(level)->coordX_SP[j] * 4.0*PI /
            (real)gridX) + cos(para->getParH(level)->coordZ_SP[j] * 4.0*PI /
                (real)gridZ)))*(real)(gridZ) / (real)(gridX);

        para->getParH(level)->vy_SP[j] = (real)0.0;
        para->getParH(level)->vx_SP[j] =
            (real)(para->getVelocity()*sin((para->getParH(level)->coordX_SP[j] * 2.0*PI / (real)gridX))*cos(para->getParH(level)->coordZ_SP[j] * 2.0*PI / (real)gridZ)) + uAdvect * (1.0 + para->getParH(level)->rho_SP[j]);
        para->getParH(level)->vz_SP[j] =
            (real)(-para->getVelocity()*cos((para->getParH(level)->coordX_SP[j] * 2.0*PI / (real)gridX))*sin(para->getParH(level)->coordZ_SP[j] * 2.0*PI / (real)gridZ));// *(real)(gridZ) / (real)(gridX);

       //para->getParH(level)->vx_SP[j] = para->getVelocity();//0.0f;//0.035f;
       //para->getParH(level)->vy_SP[j] = 0.0f;//para->getVelocity();//0.0f;
       //para->getParH(level)->vz_SP[j] = 0.0f;
       //para->getParH(level)->rho_SP[j] = 0.0f;
       //para->getParH(level)->press_SP[j] = 0.0f;
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
