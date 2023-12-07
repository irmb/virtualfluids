//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//=======================================================================================
#include "GridProvider.h"

#include "GridReaderFiles/GridReader.h"
#include "GridReaderGenerator/GridGenerator.h"
#include <Parameter/Parameter.h>

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include "Cuda/CudaMemoryManager.h"

std::shared_ptr<GridProvider> GridProvider::makeGridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager, vf::parallel::Communicator& communicator)
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

        if (para->getCalcMean()) {
            para->getParH(level)->meanVelocityInXdirection[pos] = 0.0f;
            para->getParH(level)->meanVelocityInYdirection[pos] = 0.0f;
            para->getParH(level)->meanVelocityInZdirection[pos] = 0.0f;
            para->getParH(level)->meanDensity[pos] = 0.0f;
            para->getParH(level)->meanPressure[pos] = 0.0f;
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
}


void GridProvider::setVelocitySizePerLevel(int level, int sizePerLevel) const
{
    para->getParH(level)->velocityBC.numberOfBCnodes = sizePerLevel;
    para->getParD(level)->velocityBC.numberOfBCnodes = sizePerLevel;
}

void GridProvider::setOutflowSizePerLevel(int level, int sizePerLevel) const
{
    para->getParH(level)->outflowBC.numberOfBCnodes = sizePerLevel;
    para->getParD(level)->outflowBC.numberOfBCnodes = sizePerLevel;
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
