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
//! \file ScalingUtilities.h
//! \ingroup LBM/GPUHelperFunctions
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================
#ifndef SCALING_HELPER_FUNCTIONS_H
#define SCALING_HELPER_FUNCTIONS_H

#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

namespace vf::gpu
{

__device__ __inline__ void calculateMomentsOnSourceNodes(Distributions27 &dist, real &omega, unsigned int &k_000,
                                                         unsigned int &k_M00, unsigned int &k_0M0, unsigned int &k_00M,
                                                         unsigned int &k_MM0, unsigned int &k_M0M, unsigned int &k_0MM,
                                                         unsigned int &k_MMM, real &drho, real &velocityX,
                                                         real &velocityY, real &velocityZ, real &kxyFromfcNEQ,
                                                         real &kyzFromfcNEQ, real &kxzFromfcNEQ, real &kxxMyyFromfcNEQ,
                                                         real &kxxMzzFromfcNEQ)
{
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set local distributions (f's) on source nodes:
    //!
    real f_000 = (dist.f[DIR_000])[k_000];
    real f_P00 = (dist.f[DIR_P00])[k_000];
    real f_M00 = (dist.f[DIR_M00])[k_M00];
    real f_0P0 = (dist.f[DIR_0P0])[k_000];
    real f_0M0 = (dist.f[DIR_0M0])[k_0M0];
    real f_00P = (dist.f[DIR_00P])[k_000];
    real f_00M = (dist.f[DIR_00M])[k_00M];
    real f_PP0 = (dist.f[DIR_PP0])[k_000];
    real f_MM0 = (dist.f[DIR_MM0])[k_MM0];
    real f_PM0 = (dist.f[DIR_PM0])[k_0M0];
    real f_MP0 = (dist.f[DIR_MP0])[k_M00];
    real f_P0P = (dist.f[DIR_P0P])[k_000];
    real f_M0M = (dist.f[DIR_M0M])[k_M0M];
    real f_P0M = (dist.f[DIR_P0M])[k_00M];
    real f_M0P = (dist.f[DIR_M0P])[k_M00];
    real f_0PP = (dist.f[DIR_0PP])[k_000];
    real f_0MM = (dist.f[DIR_0MM])[k_0MM];
    real f_0PM = (dist.f[DIR_0PM])[k_00M];
    real f_0MP = (dist.f[DIR_0MP])[k_0M0];
    real f_PPP = (dist.f[DIR_PPP])[k_000];
    real f_MPP = (dist.f[DIR_MPP])[k_M00];
    real f_PMP = (dist.f[DIR_PMP])[k_0M0];
    real f_MMP = (dist.f[DIR_MMP])[k_MM0];
    real f_PPM = (dist.f[DIR_PPM])[k_00M];
    real f_MPM = (dist.f[DIR_MPM])[k_M0M];
    real f_PMM = (dist.f[DIR_PMM])[k_0MM];
    real f_MMM = (dist.f[DIR_MMM])[k_MMM];

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    drho = ((((f_PPP + f_MMM) + (f_MPM + f_PMP)) + ((f_MPP + f_PMM) + (f_MMP + f_PPM))) +
            (((f_0MP + f_0PM) + (f_0MM + f_0PP)) + ((f_M0P + f_P0M) + (f_M0M + f_P0P)) +
             ((f_MP0 + f_PM0) + (f_MM0 + f_PP0))) +
            ((f_M00 + f_P00) + (f_0M0 + f_0P0) + (f_00M + f_00P))) +
           f_000;

    real oneOverRho = c1o1 / (c1o1 + drho);

    velocityX = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_PMM - f_MPP) + (f_PPM - f_MMP))) +
                 (((f_P0M - f_M0P) + (f_P0P - f_M0M)) + ((f_PM0 - f_MP0) + (f_PP0 - f_MM0))) + (f_P00 - f_M00)) *
                oneOverRho;
    velocityY = ((((f_PPP - f_MMM) + (f_MPM - f_PMP)) + ((f_MPP - f_PMM) + (f_PPM - f_MMP))) +
                 (((f_0PM - f_0MP) + (f_0PP - f_0MM)) + ((f_MP0 - f_PM0) + (f_PP0 - f_MM0))) + (f_0P0 - f_0M0)) *
                oneOverRho;
    velocityZ = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_MPP - f_PMM) + (f_MMP - f_PPM))) +
                 (((f_0MP - f_0PM) + (f_0PP - f_0MM)) + ((f_M0P - f_P0M) + (f_P0P - f_M0M))) + (f_00P - f_00M)) *
                oneOverRho;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Calculate second order moments for interpolation
    //!
    // example: kxxMzz: moment, second derivative in x direction minus the second derivative in z direction
    kxyFromfcNEQ = -c3o1 * omega *
                   ((f_MM0 + f_MMM + f_MMP - f_MP0 - f_MPM - f_MPP - f_PM0 - f_PMM - f_PMP + f_PP0 + f_PPM + f_PPP) /
                    (c1o1 + drho) -
                    ((velocityX * velocityY)));
    kyzFromfcNEQ = -c3o1 * omega *
                   ((f_0MM + f_PMM + f_MMM - f_0MP - f_PMP - f_MMP - f_0PM - f_PPM - f_MPM + f_0PP + f_PPP + f_MPP) /
                    (c1o1 + drho) -
                    ((velocityY * velocityZ)));
    kxzFromfcNEQ = -c3o1 * omega *
                   ((f_M0M + f_MMM + f_MPM - f_M0P - f_MMP - f_MPP - f_P0M - f_PMM - f_PPM + f_P0P + f_PMP + f_PPP) /
                    (c1o1 + drho) -
                    ((velocityX * velocityZ)));
    kxxMyyFromfcNEQ = -c3o2 * omega *
                      ((f_M0M + f_M00 + f_M0P - f_0MM - f_0M0 - f_0MP - f_0PM - f_0P0 - f_0PP + f_P0M + f_P00 + f_P0P) /
                       (c1o1 + drho) -
                       ((velocityX * velocityX - velocityY * velocityY)));
    kxxMzzFromfcNEQ = -c3o2 * omega *
                      ((f_MM0 + f_M00 + f_MP0 - f_0MM - f_0MP - f_00M - f_00P - f_0PM - f_0PP + f_PM0 + f_P00 + f_PP0) /
                       (c1o1 + drho) -
                       ((velocityX * velocityX - velocityZ * velocityZ)));
}

} // namespace vf::gpu

#endif
