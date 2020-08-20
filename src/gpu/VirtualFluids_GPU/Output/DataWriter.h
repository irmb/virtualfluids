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
//! \file DataWriter.h
//! \ingroup Output
//! \author Martin Schoenherr
//=======================================================================================
#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

//! \brief Class forwarding for CudaMemoryManager and Parameter
class Parameter;
class CudaMemoryManager;

//! \class FileWriter
//! \brief manages the VTK output
class DataWriter
{
public:
	//! Class default constructor
	VIRTUALFLUIDS_GPU_EXPORT DataWriter() {}
	//! Class destructor
	virtual VIRTUALFLUIDS_GPU_EXPORT ~DataWriter() {}

	//! \brief write the initialization step to VTK file(s)
	//! \param para instance of classParameter
	//! \param cudaManager instance of class CudaMemoryManager
	virtual void VIRTUALFLUIDS_GPU_EXPORT writeInit(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager) = 0;
	//! \brief write time step to VTK file(s)
	//! \param para instance of classParameter
	//! \param timestep of the simulation
	virtual void VIRTUALFLUIDS_GPU_EXPORT writeTimestep(SPtr<Parameter> para, uint timestep) = 0;
};
#endif
