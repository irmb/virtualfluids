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
//! \file CudaKernelManager.h
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
#ifndef CudaKernelManager_H
#define CudaKernelManager_H

#include <memory>
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

//! \brief Class forwarding for Parameter
class Parameter;

//! \class CudaKernelManager
//! \brief manage the cuda kernel calls
class VIRTUALFLUIDS_GPU_EXPORT CudaKernelManager
{
public:
	//! \brief makes an object of CudaKernelManager
	//! \param para shared pointer to instance of class Parameter
    static SPtr<CudaKernelManager> make(std::shared_ptr<Parameter> parameter);
    
	//! \brief calls the device function of the lattice Boltzmann kernel
	void runLBMKernel(int level);

	//! \brief calls the device function of the velocity boundary condition
    void runVelocityBCKernel(int level);

	//! \brief calls the device function of the geometry boundary condition
	void runGeoBCKernelPost(int level);

	//! \brief calls the device function of the slip boundary condition
	void runSlipBCKernel(int level);

	//! \brief calls the device function of the no-slip boundary condition
	void runNoSlipBCKernel(int level);

	//! \brief calls the device function of the pressure boundary condition
	void runPressureBCKernel(int level);

	//! \brief calls the device function of the stress wall model
	void runStressWallModelKernel(int level);


    //! \brief calls the device function that calculates the macroscopic values
    void calculateMacroscopicValues(int level);


private:
	//! Class constructor
	//! \param parameter shared pointer to instance of class Parameter
	CudaKernelManager(SPtr<Parameter> parameter);
	//! Class copy constructor
	//! \param CudaKernelManager is a reference to CudaKernelManager object
	CudaKernelManager(const CudaKernelManager&);

	//! \property para is a shared pointer to an object of Parameter
	SPtr<Parameter> para;

};
#endif
