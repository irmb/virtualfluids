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
//! \file Simulation.h
//! \ingroup LBM
//! \author Martin Schoenherr
//=======================================================================================
#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <PointerDefinitions.h>
#include <VirtualFluidsDefinitions.h>
#include "Output/LogWriter.hpp"
#include "VirtualFluids_GPU_export.h"

//! \brief Class forwarding for CudaMemoryManager, Parameter, GridProvider and DataWriter
class CudaMemoryManager;
class Parameter;
class GridProvider;
class DataWriter;

//! \class Simulation
//! \brief manages the preprocessing, simulations and post-processing
class VIRTUALFLUIDS_GPU_EXPORT Simulation
{
public:
	//! Class default constructor
	Simulation();
	//! Class destructor
	~Simulation();
	//! \brief includes the time loop over all LB-timesteps 
	void run();
	//! \brief initialize the lattice (incl. distribution functions)
	void init(SPtr<Parameter> para, SPtr<GridProvider> gridProvider, SPtr<DataWriter> dataWriter, SPtr<CudaMemoryManager> cudaManager);
	//! \brief frees the pinned host memory
	void free();

protected:
	//! \property output is an object of LogWriter
	LogWriter output;

	//! \brief shared pointer to parameter, gridProvider, dataWriter and cudaManager objects
	//! \property para is a shared pointer to an object of Parameter
	SPtr<Parameter> para;
	//! \property gridProvider is a shared pointer to an object of GridProvider
	SPtr<GridProvider> gridProvider;
	//! \property dataWriter is a shared pointer to an object of DataWriter
	SPtr<DataWriter> dataWriter;
	//! \property cudaManager is a shared pointer to an object of CudaMemoryManager
	SPtr<CudaMemoryManager> cudaManager;
 };
#endif
