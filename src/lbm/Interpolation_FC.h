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
#ifndef LBM_INTERPOLATION_FC_H
#define LBM_INTERPOLATION_FC_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif


#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"

#include "lbm/KernelParameter.h"
#include "lbm/Chimera.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

__host__ __device__ __inline__ void interpolate_fc(real* const f,
    real eps_new,
    real omegaC,
    real xoff,
    real yoff,
    real zoff,
    real xoff_sq,
    real yoff_sq,
    real zoff_sq,
    real drho_PPP, real vx1_PPP, real vx2_PPP, real vx3_PPP,
    real drho_MPP, real vx1_MPP, real vx2_MPP, real vx3_MPP,
    real drho_PMP, real vx1_PMP, real vx2_PMP, real vx3_PMP,
    real drho_MMP, real vx1_MMP, real vx2_MMP, real vx3_MMP,
    real drho_PPM, real vx1_PPM, real vx2_PPM, real vx3_PPM,
    real drho_MPM, real vx1_MPM, real vx2_MPM, real vx3_MPM,
    real drho_PMM, real vx1_PMM, real vx2_PMM, real vx3_PMM,
    real drho_MMM, real vx1_MMM, real vx2_MMM, real vx3_MMM,
    real kxyFromfcNEQ_PPP, real kyzFromfcNEQ_PPP, real kxzFromfcNEQ_PPP, real kxxMyyFromfcNEQ_PPP, real kxxMzzFromfcNEQ_PPP,
    real kxyFromfcNEQ_MPP, real kyzFromfcNEQ_MPP, real kxzFromfcNEQ_MPP, real kxxMyyFromfcNEQ_MPP, real kxxMzzFromfcNEQ_MPP,
    real kxyFromfcNEQ_PMP, real kyzFromfcNEQ_PMP, real kxzFromfcNEQ_PMP, real kxxMyyFromfcNEQ_PMP, real kxxMzzFromfcNEQ_PMP,
    real kxyFromfcNEQ_MMP, real kyzFromfcNEQ_MMP, real kxzFromfcNEQ_MMP, real kxxMyyFromfcNEQ_MMP, real kxxMzzFromfcNEQ_MMP,
    real kxyFromfcNEQ_PPM, real kyzFromfcNEQ_PPM, real kxzFromfcNEQ_PPM, real kxxMyyFromfcNEQ_PPM, real kxxMzzFromfcNEQ_PPM,
    real kxyFromfcNEQ_MPM, real kyzFromfcNEQ_MPM, real kxzFromfcNEQ_MPM, real kxxMyyFromfcNEQ_MPM, real kxxMzzFromfcNEQ_MPM,
    real kxyFromfcNEQ_PMM, real kyzFromfcNEQ_PMM, real kxzFromfcNEQ_PMM, real kxxMyyFromfcNEQ_PMM, real kxxMzzFromfcNEQ_PMM,
    real kxyFromfcNEQ_MMM, real kyzFromfcNEQ_MMM, real kxzFromfcNEQ_MMM, real kxxMyyFromfcNEQ_MMM, real kxxMzzFromfcNEQ_MMM)
{
    //////////////////////////////////////////////////////////////////////////
    //! - Calculate coefficients for polynomial interpolation
    //!
    // example: a_110: derivation in x and y direction
    real a_000, a_100, a_010, a_001, a_200, a_020, a_002, a_110, a_101, a_011;
    real b_000, b_100, b_010, b_001, b_200, b_020, b_002, b_110, b_101, b_011;
    real c_000, c_100, c_010, c_001, c_200, c_020, c_002, c_110, c_101, c_011;
    real d_000, d_100, d_010, d_001, d_110, d_101, d_011;

    a_000 = c1o64 * (
            c2o1 * (
            ((kxyFromfcNEQ_MMM - kxyFromfcNEQ_PPP) + (kxyFromfcNEQ_MMP - kxyFromfcNEQ_PPM)) + ((kxyFromfcNEQ_PMM - kxyFromfcNEQ_MPP) + (kxyFromfcNEQ_PMP - kxyFromfcNEQ_MPM)) + 
            ((kxzFromfcNEQ_MMM - kxzFromfcNEQ_PPP) + (kxzFromfcNEQ_PPM - kxzFromfcNEQ_MMP)) + ((kxzFromfcNEQ_PMM - kxzFromfcNEQ_MPP) + (kxzFromfcNEQ_MPM - kxzFromfcNEQ_PMP)) + 
            ((vx2_PPP + vx2_MMM) + (vx2_PPM + vx2_MMP)) - ((vx2_MPP + vx2_PMM) + (vx2_MPM + vx2_PMP)) + 
            ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MMP)) + ((vx3_PMP + vx3_MPM) - (vx3_MPP + vx3_PMM))) + 
            c8o1 * (((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) + ((vx1_MPP + vx1_PMM) + (vx1_PMP + vx1_MPM))) +
            ((kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_PPP) + (kxxMyyFromfcNEQ_MMP - kxxMyyFromfcNEQ_PPM)) + 
            ((kxxMyyFromfcNEQ_MPP - kxxMyyFromfcNEQ_PMM) + (kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_PMP)) +
            ((kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_PPP) + (kxxMzzFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM)) + 
            ((kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM) + (kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_PMP)));
    b_000 = c1o64 * (
            c2o1 * (
            ((kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MMM) + (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_MMP)) + 
            ((kxxMyyFromfcNEQ_MPP - kxxMyyFromfcNEQ_PMM) + (kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_PMP)) + 
            ((kxyFromfcNEQ_MMM - kxyFromfcNEQ_PPP) + (kxyFromfcNEQ_MMP - kxyFromfcNEQ_PPM)) + 
            ((kxyFromfcNEQ_MPP - kxyFromfcNEQ_PMM) + (kxyFromfcNEQ_MPM - kxyFromfcNEQ_PMP)) + 
            ((kyzFromfcNEQ_MMM - kyzFromfcNEQ_PPP) + (kyzFromfcNEQ_PPM - kyzFromfcNEQ_MMP)) + 
            ((kyzFromfcNEQ_PMM - kyzFromfcNEQ_MPP) + (kyzFromfcNEQ_MPM - kyzFromfcNEQ_PMP)) + 
            ((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) - ((vx1_MPM + vx1_MPP) + (vx1_PMM + vx1_PMP)) + 
            ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MMP)) + ((vx3_MPP + vx3_PMM) - (vx3_MPM + vx3_PMP))) + 
            c8o1 * (((vx2_PPP + vx2_MMM) + (vx2_PPM + vx2_MMP)) + ((vx2_MPP + vx2_PMM) + (vx2_MPM + vx2_PMP))) + 
            ((kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_PPP) + (kxxMzzFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM)) +
            ((kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_MPP) + (kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MPM)));
    c_000 = c1o64 * ( 
            c2o1 * (
            ((kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MMM) + (kxxMzzFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM)) + 
            ((kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM) + (kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MPM)) + 
            ((kxzFromfcNEQ_MMM - kxzFromfcNEQ_PPP) + (kxzFromfcNEQ_MMP - kxzFromfcNEQ_PPM)) + 
            ((kxzFromfcNEQ_MPP - kxzFromfcNEQ_PMM) + (kxzFromfcNEQ_MPM - kxzFromfcNEQ_PMP)) + 
            ((kyzFromfcNEQ_MMM - kyzFromfcNEQ_PPP) + (kyzFromfcNEQ_MMP - kyzFromfcNEQ_PPM)) + 
            ((kyzFromfcNEQ_PMM - kyzFromfcNEQ_MPP) + (kyzFromfcNEQ_PMP - kyzFromfcNEQ_MPM)) + 
            ((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_MPM + vx1_PMP) - (vx1_MPP + vx1_PMM)) + 
            ((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_MPP + vx2_PMM) - (vx2_MPM + vx2_PMP))) + 
            c8o1 * (((vx3_PPP + vx3_MMM) + (vx3_PPM + vx3_MMP)) + ((vx3_PMM + vx3_MPP) + (vx3_PMP + vx3_MPM))) +
            ((kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_PPP) + (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_MMP)) + 
            ((kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_MPP) + (kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_PMP)));

    a_100 = c1o4 * (((vx1_PPP - vx1_MMM) + (vx1_PPM - vx1_MMP)) + ((vx1_PMM - vx1_MPP) + (vx1_PMP - vx1_MPM)));
    b_100 = c1o4 * (((vx2_PPP - vx2_MMM) + (vx2_PPM - vx2_MMP)) + ((vx2_PMM - vx2_MPP) + (vx2_PMP - vx2_MPM)));
    c_100 = c1o4 * (((vx3_PPP - vx3_MMM) + (vx3_PPM - vx3_MMP)) + ((vx3_PMM - vx3_MPP) + (vx3_PMP - vx3_MPM)));

    a_200 = c1o16 * ( 
            c2o1 * (
            ((vx2_PPP + vx2_MMM) + (vx2_PPM - vx2_MPP)) + ((vx2_MMP - vx2_PMM) - (vx2_MPM + vx2_PMP)) + 
            ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MPP)) + ((vx3_MPM + vx3_PMP) - (vx3_MMP + vx3_PMM))) + 
            ((kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MMM) + (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_MMP)) + 
            ((kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_MPP) + (kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MPM)) + 
            ((kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MMM) + (kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_MMP)) + 
            ((kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_MPP) + (kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MPM)));
    b_200 = c1o8 * (
            c2o1 * (
            -((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) + ((vx1_MPP + vx1_PMM) + (vx1_MPM + vx1_PMP))) +
            ((kxyFromfcNEQ_PPP - kxyFromfcNEQ_MMM) + (kxyFromfcNEQ_PPM - kxyFromfcNEQ_MMP)) + 
            ((kxyFromfcNEQ_PMM - kxyFromfcNEQ_MPP) + (kxyFromfcNEQ_PMP - kxyFromfcNEQ_MPM)));
    c_200 = c1o8 * (
            c2o1 * (
            ((vx1_PPM + vx1_MMP) - (vx1_PPP + vx1_MMM)) + ((vx1_MPP + vx1_PMM) - (vx1_MPM + vx1_PMP))) +
            ((kxzFromfcNEQ_PPP - kxzFromfcNEQ_MMM) + (kxzFromfcNEQ_PPM - kxzFromfcNEQ_MMP)) + 
            ((kxzFromfcNEQ_PMM - kxzFromfcNEQ_MPP) + (kxzFromfcNEQ_PMP - kxzFromfcNEQ_MPM)));

    a_010 = c1o4 * (((vx1_PPP - vx1_MMM) + (vx1_PPM - vx1_MMP)) + ((vx1_MPP - vx1_PMM) + (vx1_MPM - vx1_PMP)));
    b_010 = c1o4 * (((vx2_PPP - vx2_MMM) + (vx2_PPM - vx2_MMP)) + ((vx2_MPP - vx2_PMM) + (vx2_MPM - vx2_PMP)));
    c_010 = c1o4 * (((vx3_PPP - vx3_MMM) + (vx3_PPM - vx3_MMP)) + ((vx3_MPP - vx3_PMM) + (vx3_MPM - vx3_PMP)));

    a_020 = c1o8 * (
            c2o1 * (-((vx2_PPP + vx2_MMM) + (vx2_MMP + vx2_PPM)) + ((vx2_MPP + vx2_PMM) + (vx2_MPM + vx2_PMP))) +
            ((kxyFromfcNEQ_PPP - kxyFromfcNEQ_MMM) + (kxyFromfcNEQ_PPM - kxyFromfcNEQ_MMP)) + 
            ((kxyFromfcNEQ_MPP - kxyFromfcNEQ_PMM) + (kxyFromfcNEQ_MPM - kxyFromfcNEQ_PMP)));
    b_020 = c1o16 * (
            c2o1 * (
            ((kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_PPP) + (kxxMyyFromfcNEQ_MMP - kxxMyyFromfcNEQ_PPM)) +
            ((kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_MPP) + (kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MPM)) +
            ((vx1_PPP + vx1_MMM) + (vx1_PPM + vx1_MMP)) - ((vx1_MPP + vx1_PMM) + (vx1_PMP + vx1_MPM)) + 
            ((vx3_PPP + vx3_MMM) - (vx3_PPM + vx3_MMP)) + ((vx3_MPP + vx3_PMM) - (vx3_MPM + vx3_PMP))) +
            ((kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MMM) + (kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_MMP)) + 
            ((kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM) + (kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_PMP)));
    c_020 = c1o8 * (
            c2o1 * (((vx2_MMP + vx2_PPM) - (vx2_PPP + vx2_MMM)) + ((vx2_PMP + vx2_MPM) - (vx2_MPP + vx2_PMM))) +
            ((kyzFromfcNEQ_PPP - kyzFromfcNEQ_MMM) + (kyzFromfcNEQ_PPM - kyzFromfcNEQ_MMP)) +
            ((kyzFromfcNEQ_MPP - kyzFromfcNEQ_PMM) + (kyzFromfcNEQ_MPM - kyzFromfcNEQ_PMP)));

    a_001 = c1o4 * (((vx1_PPP - vx1_MMM) + (vx1_MMP - vx1_PPM)) + ((vx1_MPP - vx1_PMM) + (vx1_PMP - vx1_MPM)));
    b_001 = c1o4 * (((vx2_PPP - vx2_MMM) + (vx2_MMP - vx2_PPM)) + ((vx2_MPP - vx2_PMM) + (vx2_PMP - vx2_MPM)));
    c_001 = c1o4 * (((vx3_PPP - vx3_MMM) + (vx3_MMP - vx3_PPM)) + ((vx3_MPP - vx3_PMM) + (vx3_PMP - vx3_MPM)));

    a_002 = c1o8 * (
            c2o1 * (((vx3_PPM + vx3_MMP) - (vx3_PPP + vx3_MMM)) + ((vx3_MPP + vx3_PMM) - (vx3_PMP + vx3_MPM))) +
                    ((kxzFromfcNEQ_PPP - kxzFromfcNEQ_MMM) + (kxzFromfcNEQ_MMP - kxzFromfcNEQ_PPM)) +
                    ((kxzFromfcNEQ_PMP - kxzFromfcNEQ_MPM) + (kxzFromfcNEQ_MPP - kxzFromfcNEQ_PMM)));
    b_002 = c1o8 * (
            c2o1 * (((vx3_PPM + vx3_MMP) - (vx3_PPP + vx3_MMM)) + ((vx3_MPM + vx3_PMP) - (vx3_PMM + vx3_MPP))) + 
                    ((kyzFromfcNEQ_PPP - kyzFromfcNEQ_MMM) + (kyzFromfcNEQ_MMP - kyzFromfcNEQ_PPM)) + 
                    ((kyzFromfcNEQ_PMP - kyzFromfcNEQ_MPM) + (kyzFromfcNEQ_MPP - kyzFromfcNEQ_PMM)));
    c_002 = c1o16 * (
            c2o1 * (
            ((kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_PPP) + (kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_MMP)) + 
            ((kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_PMP) + (kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_MPP)) + 
            ((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_MPM + vx1_PMP) - (vx1_PMM + vx1_MPP)) + 
            ((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_PMM + vx2_MPP) - (vx2_MPM + vx2_PMP))) + 
            ((kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MMM) + (kxxMyyFromfcNEQ_MMP - kxxMyyFromfcNEQ_PPM)) +
            ((kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MPM) + (kxxMyyFromfcNEQ_MPP - kxxMyyFromfcNEQ_PMM)));

    a_110 = c1o2 * (((vx1_PPP + vx1_MMM) + (vx1_MMP + vx1_PPM)) - ((vx1_MPM + vx1_PMP) + (vx1_PMM + vx1_MPP)));
    b_110 = c1o2 * (((vx2_PPP + vx2_MMM) + (vx2_MMP + vx2_PPM)) - ((vx2_MPM + vx2_PMP) + (vx2_PMM + vx2_MPP)));
    c_110 = c1o2 * (((vx3_PPP + vx3_MMM) + (vx3_MMP + vx3_PPM)) - ((vx3_MPM + vx3_PMP) + (vx3_PMM + vx3_MPP)));

    a_101 = c1o2 * (((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_MPM + vx1_PMP) - (vx1_PMM + vx1_MPP)));
    b_101 = c1o2 * (((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_MPM + vx2_PMP) - (vx2_PMM + vx2_MPP)));
    c_101 = c1o2 * (((vx3_PPP + vx3_MMM) - (vx3_MMP + vx3_PPM)) + ((vx3_MPM + vx3_PMP) - (vx3_PMM + vx3_MPP)));
    
    a_011 = c1o2 * (((vx1_PPP + vx1_MMM) - (vx1_MMP + vx1_PPM)) + ((vx1_PMM + vx1_MPP) - (vx1_MPM + vx1_PMP)));
    b_011 = c1o2 * (((vx2_PPP + vx2_MMM) - (vx2_MMP + vx2_PPM)) + ((vx2_PMM + vx2_MPP) - (vx2_MPM + vx2_PMP)));
    c_011 = c1o2 * (((vx3_PPP + vx3_MMM) - (vx3_MMP + vx3_PPM)) + ((vx3_PMM + vx3_MPP) - (vx3_MPM + vx3_PMP)));

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    real kxyAverage    = c0o1;
    real kyzAverage    = c0o1;
    real kxzAverage    = c0o1;
    real kxxMyyAverage = c0o1;
    real kxxMzzAverage = c0o1;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!- Calculate coefficients for the polynomial interpolation of the pressure
    //! 
    real LaplaceRho = 
        ((xoff != c0o1) || (yoff != c0o1) || (zoff != c0o1))
        ? c0o1 : -c3o1 * (a_100 * a_100 + b_010 * b_010 + c_001 * c_001) - c6o1 * (b_100 * a_010 + c_100 * a_001 + c_010 * b_001);
    d_000 =  c1o8 * ((((drho_PPP + drho_MMM) + (drho_PPM + drho_MMP)) + ((drho_PMM + drho_MPP) + (drho_PMP + drho_MPM))) - c2o1 * LaplaceRho);
    d_100 = c1o4 * (((drho_PPP - drho_MMM) + (drho_PPM - drho_MMP)) + ((drho_PMM - drho_MPP) + (drho_PMP - drho_MPM)));
    d_010 = c1o4 * (((drho_PPP - drho_MMM) + (drho_PPM - drho_MMP)) + ((drho_MPP - drho_PMM) + (drho_MPM - drho_PMP)));
    d_001 = c1o4 * (((drho_PPP - drho_MMM) + (drho_MMP - drho_PPM)) + ((drho_MPP - drho_PMM) + (drho_PMP - drho_MPM)));
    d_110 = c1o2 * (((drho_PPP + drho_MMM) + (drho_PPM + drho_MMP)) - ((drho_PMM + drho_MPP) + (drho_PMP + drho_MPM)));
    d_101 = c1o2 * (((drho_PPP + drho_MMM) - (drho_PPM + drho_MMP)) + ((drho_PMP + drho_MPM) - (drho_PMM + drho_MPP)));
    d_011 = c1o2 * (((drho_PPP + drho_MMM) - (drho_PPM + drho_MMP)) + ((drho_PMM + drho_MPP) - (drho_PMP + drho_MPM)));


    //////////////////////////////////////////////////////////////////////////
    //! - Extrapolation for refinement in to the wall (polynomial coefficients)
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // x------x
    // |      |
    // |   ---+--->X
    // |      |  |
    // x------x  |
    //          offset-vector
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    a_000 = a_000 + xoff * a_100 + yoff * a_010 + zoff * a_001 + xoff_sq * a_200 + yoff_sq * a_020 + zoff_sq * a_002 +
            xoff * yoff * a_110 + xoff * zoff * a_101 + yoff * zoff * a_011;
    a_100 = a_100 + c2o1 * xoff * a_200 + yoff * a_110 + zoff * a_101;
    a_010 = a_010 + c2o1 * yoff * a_020 + xoff * a_110 + zoff * a_011;
    a_001 = a_001 + c2o1 * zoff * a_002 + xoff * a_101 + yoff * a_011;
    b_000 = b_000 + xoff * b_100 + yoff * b_010 + zoff * b_001 + xoff_sq * b_200 + yoff_sq * b_020 + zoff_sq * b_002 +
            xoff * yoff * b_110 + xoff * zoff * b_101 + yoff * zoff * b_011;
    b_100 = b_100 + c2o1 * xoff * b_200 + yoff * b_110 + zoff * b_101;
    b_010 = b_010 + c2o1 * yoff * b_020 + xoff * b_110 + zoff * b_011;
    b_001 = b_001 + c2o1 * zoff * b_002 + xoff * b_101 + yoff * b_011;
    c_000 = c_000 + xoff * c_100 + yoff * c_010 + zoff * c_001 + xoff_sq * c_200 + yoff_sq * c_020 + zoff_sq * c_002 +
            xoff * yoff * c_110 + xoff * zoff * c_101 + yoff * zoff * c_011;
    c_100 = c_100 + c2o1 * xoff * c_200 + yoff * c_110 + zoff * c_101;
    c_010 = c_010 + c2o1 * yoff * c_020 + xoff * c_110 + zoff * c_011;
    c_001 = c_001 + c2o1 * zoff * c_002 + xoff * c_101 + yoff * c_011;
    d_000 = d_000 + xoff * d_100 + yoff * d_010 + zoff * d_001 + 
            xoff * yoff * d_110 + xoff * zoff * d_101 + yoff * zoff * d_011;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set all moments to zero
    //!
    real m_111 = c0o1;
    real m_211 = c0o1;
    real m_011 = c0o1;
    real m_121 = c0o1;
    real m_101 = c0o1;
    real m_112 = c0o1;
    real m_110 = c0o1;
    real m_221 = c0o1;
    real m_001 = c0o1;
    real m_201 = c0o1;
    real m_021 = c0o1;
    real m_212 = c0o1;
    real m_010 = c0o1;
    real m_210 = c0o1;
    real m_012 = c0o1;
    real m_122 = c0o1;
    real m_100 = c0o1;
    real m_120 = c0o1;
    real m_102 = c0o1;
    real m_222 = c0o1;
    real m_022 = c0o1;
    real m_202 = c0o1;
    real m_002 = c0o1;
    real m_220 = c0o1;
    real m_020 = c0o1;
    real m_200 = c0o1;
    real m_000 = c0o1;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Define aliases to use the same variable for the distributions (f's):
    //!
    real& f_000 = m_111;
    real& f_P00 = m_211;
    real& f_M00 = m_011;
    real& f_0P0 = m_121;
    real& f_0M0 = m_101;
    real& f_00P = m_112;
    real& f_00M = m_110;
    real& f_PP0 = m_221;
    real& f_MM0 = m_001;
    real& f_PM0 = m_201;
    real& f_MP0 = m_021;
    real& f_P0P = m_212;
    real& f_M0M = m_010;
    real& f_P0M = m_210;
    real& f_M0P = m_012;
    real& f_0PP = m_122;
    real& f_0MM = m_100;
    real& f_0PM = m_120;
    real& f_0MP = m_102;
    real& f_PPP = m_222;
    real& f_MPP = m_022;
    real& f_PMP = m_202;
    real& f_MMP = m_002;
    real& f_PPM = m_220;
    real& f_MPM = m_020;
    real& f_PMM = m_200;
    real& f_MMM = m_000;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Declare local variables for destination nodes
    //!
    real vvx, vvy, vvz, vx_sq, vy_sq, vz_sq;
    real mxxPyyPzz, mxxMyy, mxxMzz, mxxyPyzz, mxxyMyzz, mxxzPyyz, mxxzMyyz, mxyyPxzz, mxyyMxzz;
    real useNEQ = c1o1; // zero; //one;   //.... one = on ..... zero = off
    real press;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    press = d_000;
    vvx   = a_000;
    vvy   = b_000;
    vvz   = c_000;

    m_000 = press; // m_000 is press, if drho is interpolated directly

    vx_sq = vvx * vvx;
    vy_sq = vvy * vvy;
    vz_sq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    mxxPyyPzz = m_000;

    mxxMyy = -c2o3 * ((a_100 - b_010) + kxxMyyAverage) * eps_new / omegaC * (c1o1 + press);
    mxxMzz = -c2o3 * ((a_100 - c_001) + kxxMzzAverage) * eps_new / omegaC * (c1o1 + press);

    m_011 = -c1o3 * ((b_001 + c_010) + kyzAverage) * eps_new / omegaC * (c1o1 + press);
    m_101 = -c1o3 * ((a_001 + c_100) + kxzAverage) * eps_new / omegaC * (c1o1 + press);
    m_110 = -c1o3 * ((a_010 + b_100) + kxyAverage) * eps_new / omegaC * (c1o1 + press);

    m_200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m_020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m_002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    // linear combinations for third order moments
    m_111 = c0o1;

    mxxyPyzz = c0o1;
    mxxyMyzz = c0o1;
    mxxzPyyz = c0o1;
    mxxzMyyz = c0o1;
    mxyyPxzz = c0o1;
    mxyyMxzz = c0o1;

    m_210 = ( mxxyMyzz + mxxyPyzz) * c1o2;
    m_012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
    m_201 = ( mxxzMyyz + mxxzPyyz) * c1o2;
    m_021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
    m_120 = ( mxyyMxzz + mxyyPxzz) * c1o2;
    m_102 = (-mxyyMxzz + mxyyPxzz) * c1o2;

    // fourth order moments
    m_022 = m_000 * c1o9;
    m_202 = m_022;
    m_220 = m_022;

    // fifth order moments

    // sixth order moments
    m_222 = m_000 * c1o27;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardInverseChimeraWithK(m_000, m_100, m_200, vvx, vx_sq, c1o1, c1o1);
    backwardChimera(            m_010, m_110, m_210, vvx, vx_sq);
    backwardInverseChimeraWithK(m_020, m_120, m_220, vvx, vx_sq, c3o1, c1o3);
    backwardChimera(            m_001, m_101, m_201, vvx, vx_sq);
    backwardChimera(            m_011, m_111, m_211, vvx, vx_sq);
    backwardChimera(            m_021, m_121, m_221, vvx, vx_sq);
    backwardInverseChimeraWithK(m_002, m_102, m_202, vvx, vx_sq, c3o1, c1o3);
    backwardChimera(            m_012, m_112, m_212, vvx, vx_sq);
    backwardInverseChimeraWithK(m_022, m_122, m_222, vvx, vx_sq, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardInverseChimeraWithK(m_000, m_010, m_020, vvy, vy_sq, c6o1, c1o6);
    backwardChimera(            m_001, m_011, m_021, vvy, vy_sq);
    backwardInverseChimeraWithK(m_002, m_012, m_022, vvy, vy_sq, c18o1, c1o18);
    backwardInverseChimeraWithK(m_100, m_110, m_120, vvy, vy_sq, c3o2, c2o3);
    backwardChimera(            m_101, m_111, m_121, vvy, vy_sq);
    backwardInverseChimeraWithK(m_102, m_112, m_122, vvy, vy_sq, c9o2, c2o9);
    backwardInverseChimeraWithK(m_200, m_210, m_220, vvy, vy_sq, c6o1, c1o6);
    backwardChimera(            m_201, m_211, m_221, vvy, vy_sq);
    backwardInverseChimeraWithK(m_202, m_212, m_222, vvy, vy_sq, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardInverseChimeraWithK(m_000, m_001, m_002, vvz, vz_sq, c36o1, c1o36);
    backwardInverseChimeraWithK(m_010, m_011, m_012, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_020, m_021, m_022, vvz, vz_sq, c36o1, c1o36);
    backwardInverseChimeraWithK(m_100, m_101, m_102, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_110, m_111, m_112, vvz, vz_sq, c9o4,  c4o9);
    backwardInverseChimeraWithK(m_120, m_121, m_122, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_200, m_201, m_202, vvz, vz_sq, c36o1, c1o36);
    backwardInverseChimeraWithK(m_210, m_211, m_212, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_220, m_221, m_222, vvz, vz_sq, c36o1, c1o36);

    f[DIR_000] = f_000;
    f[DIR_P00] = f_P00;
    f[DIR_M00] = f_M00;
    f[DIR_0P0] = f_0P0;
    f[DIR_0M0] = f_0M0;
    f[DIR_00P] = f_00P;
    f[DIR_00M] = f_00M;
    f[DIR_PP0] = f_PP0;
    f[DIR_MM0] = f_MM0;
    f[DIR_PM0] = f_PM0;
    f[DIR_MP0] = f_MP0;
    f[DIR_P0P] = f_P0P;
    f[DIR_M0M] = f_M0M;
    f[DIR_P0M] = f_P0M;
    f[DIR_M0P] = f_M0P;
    f[DIR_0PP] = f_0PP;
    f[DIR_0MM] = f_0MM;
    f[DIR_0PM] = f_0PM;
    f[DIR_0MP] = f_0MP;
    f[DIR_PPP] = f_PPP;
    f[DIR_MPP] = f_MPP;
    f[DIR_PMP] = f_PMP;
    f[DIR_MMP] = f_MMP;
    f[DIR_PPM] = f_PPM;
    f[DIR_MPM] = f_MPM;
    f[DIR_PMM] = f_PMM;
    f[DIR_MMM] = f_MMM;
}

}

#endif
