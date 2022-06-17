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
//! \file GridScalingKernelManager.h
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GridScalingKernelManager_H
#define GridScalingKernelManager_H

#include <memory>
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

//! \brief Class forwarding for Parameter
class Parameter;
class CudaMemoryManager;

//! \class GridScalingKernelManager
//! \brief manage the cuda kernel calls
class VIRTUALFLUIDS_GPU_EXPORT GridScalingKernelManager
{
public:
	//! \brief makes an object of GridScalingKernelManager
	//! \param para shared pointer to instance of class Parameter
    static SPtr<GridScalingKernelManager> make(std::shared_ptr<Parameter> parameter);
    
	//! \brief calls the device function of the fine to coarse grid interpolation kernel
	void runFineToCoarseKernelLB(int level);

	//! \brief calls the device function of the fine to coarse grid interpolation kernel (advection diffusion)
    void runFineToCoarseKernelAD(int level);

	//! \brief calls the device function of the coarse to fine grid interpolation kernel
	void runCoarseToFineKernelLB(int level);

	//! \brief calls the device function of the coarse to fine grid interpolation kernel (advection diffusion)
    void runCoarseToFineKernelAD(int level);

private:
	//! Class constructor
	//! \param parameter shared pointer to instance of class Parameter
	GridScalingKernelManager(SPtr<Parameter> parameter);
	//! Class copy constructor
	//! \param GridScalingKernelManager is a reference to GridScalingKernelManager object
	GridScalingKernelManager(const GridScalingKernelManager&);

	//! \property para is a shared pointer to an object of Parameter
	SPtr<Parameter> para;

};
#endif
