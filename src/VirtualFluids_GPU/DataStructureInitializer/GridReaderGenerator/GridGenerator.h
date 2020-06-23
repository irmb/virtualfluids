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
//! \file GridGenerator.h
//! \ingroup DataStructureInitializer
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GridReaderGenerator_H
#define GridReaderGenerator_H

#include "../GridProvider.h"
#include <VirtualFluidsDefinitions.h>
#include "LBM/LB.h"

//! \brief Class forwarding for Parameter and GridBuilder
class Parameter;
class GridBuilder;

//! \class GridGenerator derived class of GridProvider
//! \brief mapping the grid of grid generator to data structure for simulation
class GridGenerator
	: public GridProvider
{
private:
	//! \brief string vector with channel direction
	std::vector<std::string> channelDirections;
	//! \brief string vector with channel direction
	std::vector<std::string> channelBoundaryConditions;

	//! \brief shared pointer to GridBuilder object
	//! \property builder is a shared pointer to an object of GridBuilder
	SPtr<GridBuilder> builder;

public:
	//! Class constructor
	//! \param builder shared pointer to instance of GridBuilder
	//! \param para shared pointer to instance of classParameter
	//! \param cudaManager shared pointer to instance of class CudaMemoryManager
	VF_PUBLIC GridGenerator(SPtr<GridBuilder> builder, SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager);
	//! Class default destructor
	VF_PUBLIC virtual ~GridGenerator();

	//! \brief allocates and initialized the data structures for Coordinates and node types
	void allocArrays_CoordNeighborGeo() override;
	//! \brief allocates and initialized the values at the boundary conditions
	void allocArrays_BoundaryValues() override;
	//! \brief allocates and initialized the sub-grid distances at the boundary conditions
	void allocArrays_BoundaryQs() override;
	
private:
	//! \brief verifies if there are invalid nodes, stopper nodes or wrong neighbors
	std::string verifyNeighborIndices() const;
	//! \brief verifies single neighbor index
	//! \param index type integer
	//! \param invalidNodes reference to invalid nodes
	//! \param stopperNodes reference to stopper nodes
	//! \param wrongNeighbors reference to wrong neighbors
	std::string verifyNeighborIndex(int index, int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const;
	//! \brief check the neighbors
	//! \param x,y,z lattice node position
	//! \param numberOfWrongNeihgbors reference to the number of wrong neighbors
	//! \param neighborIndex index of neighbor node
	//! \param neighborX,neighborY,neighborZ neighbor lattice node position
	//! \param direction type string
	std::string checkNeighbor(real x, real y, real z, int index, int& numberOfWrongNeihgbors, int neighborIndex, real neighborX, real neighborY, real neighborZ, std::string direction) const;
};

#endif
