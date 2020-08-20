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
//! \file CudaMemoryManager.h
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#ifndef CudamemoryManager_H
#define CudamemoryManager_H

#include <memory>
#include "Core/PointerDefinitions.h"

//! \brief Class forwarding for Parameter
class Parameter;

//! \class CudaMemoryManager
//! \brief manage the cuda memory, f.e. allocate, copy and free the memory
class VIRTUALFLUIDS_GPU_EXPORT CudaMemoryManager
{
public:
	//! \brief makes an object of CudaMemoryManager
	//! \param para shared pointer to instance of class Parameter
	static std::shared_ptr<CudaMemoryManager> make(std::shared_ptr<Parameter> parameter);
    
	//! \brief adds partial memory to set the total device memory usage
	void setMemsizeGPU(double admem, bool reset);
	//! \brief returns the actual total device memory usage
	double getMemsizeGPU();

	//! \brief allocate, copy and free the host / device memory for the coordinates
	void cudaAllocCoord();
	void cudaCopyCoord();
	void cudaFreeCoord();

	//! \brief copy macroscopic values from device to host memory
	void cudaCopyDataToHost();

	//! \brief allocate, copy and free the host / device memory for the macroscopic values
	void cudaAllocSP();
	void cudaCopySP();
	void cudaFreeSP();
    
	//! \brief allocate, copy and free the host / device memory for the neighbor in negative diagonal direction (-x,-y,-z)
	void cudaAllocNeighborWSB();
	void cudaCopyNeighborWSB();
	void cudaFreeNeighborWSB();

	//! \brief allocate, copy and free the host / device memory for the velocity boundary condition
	void cudaAllocVeloBC();
	void cudaCopyVeloBC();
	void cudaFreeVeloBC();

	//! \brief allocate, copy and free the host / device memory for the forcing
	void cudaAllocForcing();
	void cudaCopyForcingToDevice();
	void cudaCopyForcingToHost();
	void cudaFreeForcing();


private:
	//! Class constructor
	//! \param parameter shared pointer to instance of class Parameter
	CudaMemoryManager(SPtr<Parameter> parameter);
	//! Class copy constructor
	//! \param CudaMemoryManager is a reference to CudaMemoryManager object
	CudaMemoryManager(const CudaMemoryManager&);

	//! \property para is a shared pointer to an object of Parameter
	SPtr<Parameter> parameter;
	//! \property memsizeGPU stores the used device memory
	double memsizeGPU;

};
#endif
