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
//! \file GridProvider.h
//! \ingroup DataStructureInitializer
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GridReader_H
#define GridReader_H

#include <memory>

#include <VirtualFluidsDefinitions.h>
#include "Core/PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

//! \brief Class forwarding for CudaMemoryManager, GridBuilder and Parameter
class Parameter;
class GridBuilder;
class CudaMemoryManager;

//! \class GridProvider
//! \brief mapping the grid of grid generator to data structure for simulation
class VIRTUALFLUIDS_GPU_EXPORT GridProvider
{
public:
	//! \brief makes an object of GridGenerator
	//! \param builder shared pointer to instance of GridBuilder
	//! \param para shared pointer to instance of classParameter
	//! \param cudaManager shared pointer to instance of class CudaMemoryManager
	static SPtr<GridProvider> makeGridGenerator(SPtr<GridBuilder> builder, SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager);

	//! \brief allocates and initialized the data structures for Coordinates and node types
	virtual void allocArrays_CoordNeighborGeo() = 0;
	//! \brief allocates and initialized the values at the boundary conditions
	virtual void allocArrays_BoundaryValues() = 0;
	//! \brief allocates and initialized the sub-grid distances at the boundary conditions
	virtual void allocArrays_BoundaryQs() = 0;

	//! \brief allocates forces and copy them to the device
	virtual void allocAndCopyForcing();
	//! \brief clears pinned memory on host
	virtual void freeMemoryOnHost();
	//! \brief copy data from device to host
	virtual void cudaCopyDataToHost();

	//! Class default destructor
	virtual ~GridProvider() {}

protected:
	//! \brief set the number of lattice nodes and memory sizes for different data types
	//! \param numberOfNodes integer value with number of lattice nodes
	void setNumberOfNodes(const int numberOfNodes) const;
	//! \brief set the values for macroscopic values (velocities, rho) to zero
	//! \param numberOfNodes integer value with number of lattice nodes
	virtual void setInitalNodeValues(const int numberOfNodes) const;
	//! \brief set the size of velocity boundary condition
	//! \param numberOfNodes integer value with number of lattice nodes
	void setVelocitySize(int size) const;

	//! \brief shared pointer to parameter, gridProvider, dataWriter and cudaManager objects
	//! \property para is a shared pointer to an object of Parameter
	SPtr<Parameter> para;
	//! \property cudaManager is a shared pointer to an object of CudaMemoryManager
	SPtr<CudaMemoryManager> cudaMemoryManager;
};

#endif
