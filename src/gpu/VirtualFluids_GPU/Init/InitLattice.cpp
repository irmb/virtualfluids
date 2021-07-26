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
//! \file InitLattice.h
//! \ingroup Init
//! \author Martin Schoenherr
//=======================================================================================
#include "Init/InitLattice.h"

#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "PreProcessor/PreProcessor.h"
#include "Temperature/FindTemperature.h"


void initLattice(SPtr<Parameter> para, SPtr<PreProcessor> preProcessor, SPtr<CudaMemoryManager> cudaManager)
{
    for (int lev = para->getFine(); lev >= para->getCoarse(); lev--) {
        preProcessor->init(para, lev);

        CalcMacCompSP27(
            para->getParD(lev)->vx_SP, para->getParD(lev)->vy_SP, para->getParD(lev)->vz_SP, para->getParD(lev)->rho_SP,
            para->getParD(lev)->press_SP, para->getParD(lev)->geoSP, para->getParD(lev)->neighborX_SP,
            para->getParD(lev)->neighborY_SP, para->getParD(lev)->neighborZ_SP, para->getParD(lev)->size_Mat_SP,
            para->getParD(lev)->numberofthreads, para->getParD(lev)->d0SP.f[0], para->getParD(lev)->evenOrOdd);

        if (para->getCalcMedian()) {
            constexpr uint tdiff = 1;
            CalcMacMedSP27(para->getParD(lev)->vx_SP_Med, para->getParD(lev)->vy_SP_Med, para->getParD(lev)->vz_SP_Med,
                           para->getParD(lev)->rho_SP_Med, para->getParD(lev)->press_SP_Med, para->getParD(lev)->geoSP,
                           para->getParD(lev)->neighborX_SP, para->getParD(lev)->neighborY_SP,
                           para->getParD(lev)->neighborZ_SP, tdiff, para->getParD(lev)->size_Mat_SP,
                           para->getParD(lev)->numberofthreads, para->getParD(lev)->evenOrOdd);
        }
        // advection - diffusion
        if (para->getDiffOn()) {

            cudaManager->cudaAllocConc(lev);

            for (unsigned int i = 0; i < para->getParH(lev)->size_Mat_SP; i++) {
                para->getParH(lev)->Conc[i] = para->getTemperatureInit();
            }
            initTemperatur(para.get(), cudaManager.get(), lev);
        }
    }
}
