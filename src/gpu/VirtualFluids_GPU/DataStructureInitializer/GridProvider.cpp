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
//! \file GridProvider.cpp
//! \ingroup DataStructureInitializer
//! \author Martin Schoenherr
//=======================================================================================
#include "GridProvider.h"

#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"
#include "GridReaderGenerator/GridGenerator.h"
#include "GridGenerator/grid/GridBuilder/GridBuilder.h"


SPtr<GridProvider> GridProvider::makeGridGenerator(SPtr<GridBuilder> builder, SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager)
{
    return SPtr<GridProvider>(new GridGenerator(builder, para, cudaManager));
}

void GridProvider::setNumberOfNodes(const int numberOfNodes) const
{
    para->getParH()->numberOfNodes = numberOfNodes;
    para->getParD()->numberOfNodes = numberOfNodes;
    para->getParH()->mem_size_real = sizeof(real) * para->getParH()->numberOfNodes;
    para->getParH()->mem_size_int  = sizeof(uint) * para->getParH()->numberOfNodes;
    para->getParD()->mem_size_real = sizeof(real) * para->getParD()->numberOfNodes;
    para->getParD()->mem_size_int  = sizeof(uint) * para->getParD()->numberOfNodes;
}

void GridProvider::setInitalNodeValues(const int numberOfNodes) const
{
    for (int j = 1; j <= numberOfNodes; j++)
    {
        para->getParH()->rho[j] = real(0.0);
        para->getParH()->velocityY[j]  = real(0.0);
        para->getParH()->velocityX[j]  = real(0.0);
        para->getParH()->velocityZ[j]  = real(0.0);
    }
}


void GridProvider::setVelocitySize(int size) const
{
    para->getParH()->numberOfInflowBCnodes = size;
    para->getParD()->numberOfInflowBCnodes = size;
}

void GridProvider::allocAndCopyForcing()
{
    cudaMemoryManager->cudaAllocForcing();
    cudaMemoryManager->cudaCopyForcingToDevice();
}

void GridProvider::freeMemoryOnHost()
{
    cudaMemoryManager->cudaFreeCoord();
    cudaMemoryManager->cudaFreeSP();
}

void GridProvider::cudaCopyDataToHost()
{
    cudaMemoryManager->cudaCopyDataToHost();
}
