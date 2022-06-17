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
//! \file LBKernelManager.h
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

//! \class LBKernelManager
//! \brief manage the cuda kernel calls
class VIRTUALFLUIDS_GPU_EXPORT LBKernelManager
{
public:
	//! \brief makes an object of LBKernelManager
	//! \param para shared pointer to instance of class Parameter
    static SPtr<LBKernelManager> make(std::shared_ptr<Parameter> parameter);
    
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

	//! \brief calls the device function of the pressure boundary condition (pre-collision)
	void runPressureBCKernelPre(int level);

	//! \brief calls the device function of the pressure boundary condition (post-collision)
	void runPressureBCKernelPost(int level);

	//! \brief calls the device function of the outflow boundary condition
	void runOutflowBCKernel(int level);

	//! \brief calls the device function of the stress wall model
	void runStressWallModelKernel(int level);


    //! \brief calls the device function that calculates the macroscopic values
    void calculateMacroscopicValues(int level);


private:
	//! Class constructor
	//! \param parameter shared pointer to instance of class Parameter
	LBKernelManager(SPtr<Parameter> parameter);
	//! Class copy constructor
	//! \param LBKernelManager is a reference to LBKernelManager object
	LBKernelManager(const LBKernelManager&);

	//! \property para is a shared pointer to an object of Parameter
	SPtr<Parameter> para;

};
#endif
