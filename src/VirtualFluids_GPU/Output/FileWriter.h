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
//! \file FileWriter.h
//! \ingroup Output
//! \author Martin Schoenherr
//=======================================================================================
#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <vector>

#include "DataWriter.h"

//! \brief Class forwarding for CudaMemoryManager and Parameter
class Parameter;
class CudaMemoryManager;

//! \class FileWriter derived class of DataWriter
//! \brief manages the VTK output
class FileWriter : public DataWriter
{
public:
	//! Class default constructor
	VF_PUBLIC FileWriter() {}
	//! Class destructor
	VF_PUBLIC ~FileWriter() {}

	//! \brief write the initialization step to VTK file(s)
	//! \param para instance of classParameter
	//! \param cudaManager instance of class CudaMemoryManager
	void VF_PUBLIC writeInit(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager) override;
	//! \brief write time step to VTK file(s)
	//! \param para instance of classParameter
	//! \param timestep of the simulation
	void VF_PUBLIC writeTimestep(SPtr<Parameter> para, uint timestep) override;

private:
	//! \brief write binary VTK file as unstructured grid
	//! \param para instance of classParameter
	//! \param fname vector of strings with path and prefix of written files
	void VF_PUBLIC writeUnstrucuredGridLT(SPtr<Parameter> para, std::vector<std::string >& fname);
	//! \brief checks for periodic cells
	//! \param para instance of classParameter
	//! \param number 2, 1, 3 and 5 are the possible periodic neighbors
	bool VF_PUBLIC isPeriodicCell(SPtr<Parameter> para, uint number2, uint number1, uint number3, uint number5);
};
#endif