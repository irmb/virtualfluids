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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file GridGenerator.cpp
//! \ingroup DataStructureInitializer
//! \author Martin Schoenherr
//=======================================================================================
#include "GridGenerator.h"
#include <sstream>
#include <iostream>

#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Parameter/Parameter.h"
#include "GridGenerator/grid/GridBuilder/GridBuilder.h"
#include "GPU/CudaMemoryManager.h"
#include "utilities/math/Math.h"

GridGenerator::GridGenerator(SPtr<GridBuilder> builder, SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager)
{
	this->builder = builder;
    this->para = para;
    this->cudaMemoryManager = cudaManager;
}

GridGenerator::~GridGenerator()
{

}

void GridGenerator::allocArrays_CoordNeighborGeo()
{
	int numberOfNodesGlobal = 0;
	std::cout << "Number of Nodes: " << std::endl;
	
	const int numberOfNodes = builder->getNumberOfNodes(0) + 1;
	numberOfNodesGlobal += numberOfNodes;
	
	setNumberOfNodes(numberOfNodes);
	
	cudaMemoryManager->cudaAllocCoord();
    cudaMemoryManager->cudaAllocSP();
	cudaMemoryManager->cudaAllocNeighborWSB();

	builder->getNodeValues(
		para->getParH()->coordinateX,
		para->getParH()->coordinateY,
		para->getParH()->coordinateZ,
		para->getParH()->neighborX,
		para->getParH()->neighborY,
		para->getParH()->neighborZ,
		para->getParH()->neighborInverse,
		para->getParH()->typeOfGridNode,
		0);

	setInitalNodeValues(numberOfNodes);

    cudaMemoryManager->cudaCopyNeighborWSB();
    cudaMemoryManager->cudaCopySP();
    cudaMemoryManager->cudaCopyCoord();

	std::cout << "Number of Nodes: " << numberOfNodesGlobal << std::endl;
	std::cout << "-----finish Coord, Neighbor, Geo------" << std::endl;
}

void GridGenerator::allocArrays_BoundaryValues()
{
	std::cout << "------read BoundaryValues------" << std::endl;

    const auto numberOfVelocityValues = int(builder->getVelocitySize(0));
	std::cout << "size velocity: " << numberOfVelocityValues << std::endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int blocks = (numberOfVelocityValues / para->getParH()->numberofthreads) + 1;
    para->getParH()->inflowBC.kArray = blocks * para->getParH()->numberofthreads;
    para->getParD()->inflowBC.kArray = para->getParH()->inflowBC.kArray;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    para->getParH()->numberOfInflowBCnodes = numberOfVelocityValues;
    para->getParD()->numberOfInflowBCnodes = numberOfVelocityValues;

    if (numberOfVelocityValues > 1)
    {
		cudaMemoryManager->cudaAllocVeloBC();

        builder->getVelocityValues(para->getParH()->inflowBC.Vx, para->getParH()->inflowBC.Vy, para->getParH()->inflowBC.Vz, para->getParH()->inflowBC.k, 0);

		cudaMemoryManager->cudaCopyVeloBC();
    }
}


void GridGenerator::allocArrays_BoundaryQs()
{
	std::cout << "------read BoundaryQs-------" << std::endl;
    const auto numberOfVelocityNodes = int(builder->getVelocitySize(0));
    if (numberOfVelocityNodes > 0)
    {
		std::cout << "size velocity: " << numberOfVelocityNodes << std::endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //preprocessing
        real* QQ = para->getParH()->inflowBC.q27[0];
		unsigned int sizeQ = para->getParH()->numberOfInflowBCnodes;
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
        Q.q27[dirREST] = &QQ[dirREST*sizeQ];
        Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
        Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
        Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
        Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
        Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
        Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
        Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
        Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];

        builder->getVelocityQs(Q.q27, 0);

		cudaMemoryManager->cudaCopyVeloBC();
    }

	std::cout << "-----finish BoundaryQs------" << std::endl;
}




std::string GridGenerator::verifyNeighborIndices() const
{
    std::ostringstream oss;
    oss << "---------report start---------\n";
    oss << "Checking neighbor indices in grid \n";

    int invalidNodes = 0;
    int wrongNeighbors = 0;
    int stopperNodes = 0;

    for (uint index = 0; index < para->getParH()->numberOfNodes; index++)
        oss << verifyNeighborIndex(index, invalidNodes, stopperNodes, wrongNeighbors);


    oss << "invalid nodes found: " << invalidNodes << "\n";
    oss << "wrong neighbors found: " << wrongNeighbors << "\n";
    oss << "stopper nodes found : " << stopperNodes << "\n";
    oss << "---------report end---------\n";
    return oss.str();
}

std::string GridGenerator::verifyNeighborIndex(int index , int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const
{
    std::ostringstream oss;

    const int geo = para->getParH()->typeOfGridNode[index];
    if (geo == 16)
    {
        stopperNodes++;
        return "";
    }

    real x = para->getParH()->coordinateX[index];
    real y = para->getParH()->coordinateY[index];
    real z = para->getParH()->coordinateZ[index];

    real delta = para->getParH()->coordinateX[2] - para->getParH()->coordinateX[1];

    real maxX = para->getParH()->coordinateX[para->getParH()->numberOfNodes - 1] - delta;
    real maxY = para->getParH()->coordinateY[para->getParH()->numberOfNodes - 1] - delta;
    real maxZ = para->getParH()->coordinateZ[para->getParH()->numberOfNodes - 1] - delta;
    real realNeighborX = vf::Math::lessEqual(x + delta, maxX) ? x + delta : para->getParH()->coordinateX[1];
    real realNeighborY = vf::Math::lessEqual(y + delta, maxY) ? y + delta : para->getParH()->coordinateY[1];
    real realNeighborZ = vf::Math::lessEqual(z + delta, maxZ) ? z + delta : para->getParH()->coordinateZ[1];

    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborX[index], realNeighborX, y, z, "X");
    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborY[index], x, realNeighborY, z, "Y");
    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborZ[index], x, y, realNeighborZ, "Z");

    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborY[this->para->getParH()->neighborX[index]], realNeighborX, realNeighborY, z, "XY");
    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborZ[this->para->getParH()->neighborX[index]], realNeighborX, y, realNeighborZ, "XZ");
    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborZ[this->para->getParH()->neighborY[index]], x, realNeighborY, realNeighborZ, "YZ");

    oss << checkNeighbor(x, y, z, index, wrongNeighbors, this->para->getParH()->neighborZ[this->para->getParH()->neighborY[this->para->getParH()->neighborX[index]]], realNeighborX, realNeighborY, realNeighborZ, "XYZ");

    return oss.str();
}

std::string GridGenerator::checkNeighbor(real x, real y, real z, int index, int& numberOfWrongNeihgbors, int neighborIndex, real neighborX, real neighborY, real neighborZ, std::string direction) const
{
    std::ostringstream oss("");

    real neighborCoordX = para->getParH()->coordinateX[neighborIndex];
    real neighborCoordY = para->getParH()->coordinateY[neighborIndex];
    real neighborCoordZ = para->getParH()->coordinateZ[neighborIndex];

    const bool neighborValid = vf::Math::equal(neighborX, neighborCoordX) && vf::Math::equal(neighborY, neighborCoordY) && vf::Math::equal(neighborZ, neighborCoordZ);

    if (!neighborValid) {
        oss << "NeighborX invalid from: (" << x << ", " << y << ", " << z << "), index: " << index << ", "
            << direction << " neighborIndex: " << neighborIndex << 
            ", actual neighborCoords : (" << neighborCoordX << ", " << neighborCoordY << ", " << neighborCoordZ << 
            "), expected neighborCoords : (" << neighborX << ", " << neighborY << ", " << neighborZ << ")\n";
        numberOfWrongNeihgbors++;
    }
    return oss.str();
}
