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
//=======================================================================================
#ifndef LBM_SCALING_HELPER_FUNCTIONS_H
#define LBM_SCALING_HELPER_FUNCTIONS_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif


#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"

#include "lbm/KernelParameter.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

__host__ __device__ __inline__ void calculateMomentsOnSourceNodes(const real* const f, real &omega, real &drho, real &velocityX,
                                                         real &velocityY, real &velocityZ, real &kxyFromfcNEQ,
                                                         real &kyzFromfcNEQ, real &kxzFromfcNEQ, real &kxxMyyFromfcNEQ,
                                                         real &kxxMzzFromfcNEQ)
{
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set local distributions (f's) on source nodes:
    //!
    const real f_000 = f[dir::DIR_000];
    const real f_P00 = f[dir::DIR_P00];
    const real f_M00 = f[dir::DIR_M00];
    const real f_0P0 = f[dir::DIR_0P0];
    const real f_0M0 = f[dir::DIR_0M0];
    const real f_00P = f[dir::DIR_00P];
    const real f_00M = f[dir::DIR_00M];
    const real f_PP0 = f[dir::DIR_PP0];
    const real f_MM0 = f[dir::DIR_MM0];
    const real f_PM0 = f[dir::DIR_PM0];
    const real f_MP0 = f[dir::DIR_MP0];
    const real f_P0P = f[dir::DIR_P0P];
    const real f_M0M = f[dir::DIR_M0M];
    const real f_P0M = f[dir::DIR_P0M];
    const real f_M0P = f[dir::DIR_M0P];
    const real f_0PP = f[dir::DIR_0PP];
    const real f_0MM = f[dir::DIR_0MM];
    const real f_0PM = f[dir::DIR_0PM];
    const real f_0MP = f[dir::DIR_0MP];
    const real f_PPP = f[dir::DIR_PPP];
    const real f_MPP = f[dir::DIR_MPP];
    const real f_PMP = f[dir::DIR_PMP];
    const real f_MMP = f[dir::DIR_MMP];
    const real f_PPM = f[dir::DIR_PPM];
    const real f_MPM = f[dir::DIR_MPM];
    const real f_PMM = f[dir::DIR_PMM];
    const real f_MMM = f[dir::DIR_MMM];

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

    const real oneOverRho = c1o1 / (c1o1 + drho);

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

}

#endif
