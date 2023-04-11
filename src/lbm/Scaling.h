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

struct MomentsOnSourceNode {
    real drho;
    real velocityX;
    real velocityY;
    real velocityZ;
    real kxyFromfcNEQ;
    real kyzFromfcNEQ;
    real kxzFromfcNEQ;
    real kxxMyyFromfcNEQ;
    real kxxMzzFromfcNEQ;
};


__host__ __device__ __inline__ void calculateMomentsOnSourceNodes(const real* const f, const real omega, MomentsOnSourceNode &moments)
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
    moments.drho = ((((f_PPP + f_MMM) + (f_MPM + f_PMP)) + ((f_MPP + f_PMM) + (f_MMP + f_PPM))) +
            (((f_0MP + f_0PM) + (f_0MM + f_0PP)) + ((f_M0P + f_P0M) + (f_M0M + f_P0P)) +
             ((f_MP0 + f_PM0) + (f_MM0 + f_PP0))) +
            ((f_M00 + f_P00) + (f_0M0 + f_0P0) + (f_00M + f_00P))) +
           f_000;

    const real oneOverRho = c1o1 / (c1o1 + moments.drho);

    moments.velocityX = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_PMM - f_MPP) + (f_PPM - f_MMP))) +
                 (((f_P0M - f_M0P) + (f_P0P - f_M0M)) + ((f_PM0 - f_MP0) + (f_PP0 - f_MM0))) + (f_P00 - f_M00)) *
                oneOverRho;
    moments.velocityY = ((((f_PPP - f_MMM) + (f_MPM - f_PMP)) + ((f_MPP - f_PMM) + (f_PPM - f_MMP))) +
                 (((f_0PM - f_0MP) + (f_0PP - f_0MM)) + ((f_MP0 - f_PM0) + (f_PP0 - f_MM0))) + (f_0P0 - f_0M0)) *
                oneOverRho;
    moments.velocityZ = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_MPP - f_PMM) + (f_MMP - f_PPM))) +
                 (((f_0MP - f_0PM) + (f_0PP - f_0MM)) + ((f_M0P - f_P0M) + (f_P0P - f_M0M))) + (f_00P - f_00M)) *
                oneOverRho;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Calculate second order moments for interpolation
    //!
    // example: kxxMzz: moment, second derivative in x direction minus the second derivative in z direction

    moments.kxyFromfcNEQ = -c3o1 * omega *
                   ((f_MM0 + f_MMM + f_MMP - f_MP0 - f_MPM - f_MPP - f_PM0 - f_PMM - f_PMP + f_PP0 + f_PPM + f_PPP) /
                    (c1o1 + moments.drho) -
                    ((moments.velocityX * moments.velocityY)));
    moments.kyzFromfcNEQ = -c3o1 * omega *
                   ((f_0MM + f_PMM + f_MMM - f_0MP - f_PMP - f_MMP - f_0PM - f_PPM - f_MPM + f_0PP + f_PPP + f_MPP) /
                    (c1o1 + moments.drho) -
                    ((moments.velocityY * moments.velocityZ)));
    moments.kxzFromfcNEQ = -c3o1 * omega *
                   ((f_M0M + f_MMM + f_MPM - f_M0P - f_MMP - f_MPP - f_P0M - f_PMM - f_PPM + f_P0P + f_PMP + f_PPP) /
                    (c1o1 + moments.drho) -
                    ((moments.velocityX * moments.velocityZ)));
    moments.kxxMyyFromfcNEQ = -c3o2 * omega *
                      ((f_M0M + f_M00 + f_M0P - f_0MM - f_0M0 - f_0MP - f_0PM - f_0P0 - f_0PP + f_P0M + f_P00 + f_P0P) /
                       (c1o1 + moments.drho) -
                       ((moments.velocityX * moments.velocityX - moments.velocityY * moments.velocityY)));
    moments.kxxMzzFromfcNEQ = -c3o2 * omega *
                      ((f_MM0 + f_M00 + f_MP0 - f_0MM - f_0MP - f_00M - f_00P - f_0PM - f_0PP + f_PM0 + f_P00 + f_PP0) /
                       (c1o1 + moments.drho) -
                       ((moments.velocityX * moments.velocityX - moments.velocityZ * moments.velocityZ)));
}


struct Coefficients
{
    real a_000, a_100, a_010, a_001, a_200, a_020, a_002, a_110, a_101, a_011;
    real b_000, b_100, b_010, b_001, b_200, b_020, b_002, b_110, b_101, b_011;
    real c_000, c_100, c_010, c_001, c_200, c_020, c_002, c_110, c_101, c_011;
    real d_000, d_100, d_010, d_001, d_110, d_101, d_011;
    real a_111, b_111, c_111, d_111;
    real LaplaceRho;
};


__host__ __device__ __inline__ void calculateCoefficients(real xoff, real yoff, real zoff, Coefficients &coefficients, 
    const vf::lbm::MomentsOnSourceNode &moments_PPP,
    const vf::lbm::MomentsOnSourceNode &moments_MPP,
    const vf::lbm::MomentsOnSourceNode &moments_PMP,
    const vf::lbm::MomentsOnSourceNode &moments_MMP,
    const vf::lbm::MomentsOnSourceNode &moments_PPM,
    const vf::lbm::MomentsOnSourceNode &moments_MPM,
    const vf::lbm::MomentsOnSourceNode &moments_PMM,
    const vf::lbm::MomentsOnSourceNode &moments_MMM)
{
    real& a_000 = coefficients.a_000;
    real& b_000 = coefficients.b_000;
    real& c_000 = coefficients.c_000;
    real& d_000 = coefficients.d_000;

    real& a_100 = coefficients.a_100;
    real& b_100 = coefficients.b_100;
    real& c_100 = coefficients.c_100;
    real& d_100 = coefficients.d_100;

    real& a_010 = coefficients.a_010;
    real& b_010 = coefficients.b_010;
    real& c_010 = coefficients.c_010;
    real& d_010 = coefficients.d_010;

    real& a_001 = coefficients.a_001;
    real& b_001 = coefficients.b_001;
    real& c_001 = coefficients.c_001;
    real& d_001 = coefficients.d_001;

    real& d_110 = coefficients.d_110, &d_101 = coefficients.d_101, &d_011 = coefficients.d_011;
    
    real& a_200 = coefficients.a_200, &a_020 = coefficients.a_020, &a_002 = coefficients.a_002;
    real& b_200 = coefficients.b_200, &b_020 = coefficients.b_020, &b_002 = coefficients.b_002;
    real& c_200 = coefficients.c_200, &c_020 = coefficients.c_020, &c_002 = coefficients.c_002;

    real& a_110 = coefficients.a_110, &a_101 = coefficients.a_101, &a_011 = coefficients.a_011;
    real& b_110 = coefficients.b_110, &b_101 = coefficients.b_101, &b_011 = coefficients.b_011;
    real& c_110 = coefficients.c_110, &c_101 = coefficients.c_101, &c_011 = coefficients.c_011;

    real &a_111 = coefficients.a_111, &b_111 = coefficients.b_111, &c_111 = coefficients.c_111, &d_111 = coefficients.d_111;

    real &LaplaceRho = coefficients.LaplaceRho;

    const real xoff_sq = xoff * xoff;
    const real yoff_sq = yoff * yoff;
    const real zoff_sq = zoff * zoff;

    const real drho_PPP = moments_PPP.drho, vx1_PPP = moments_PPP.velocityX, vx2_PPP = moments_PPP.velocityY, vx3_PPP = moments_PPP.velocityZ;
    const real drho_MPP = moments_MPP.drho, vx1_MPP = moments_MPP.velocityX, vx2_MPP = moments_MPP.velocityY, vx3_MPP = moments_MPP.velocityZ;
    const real drho_PMP = moments_PMP.drho, vx1_PMP = moments_PMP.velocityX, vx2_PMP = moments_PMP.velocityY, vx3_PMP = moments_PMP.velocityZ;
    const real drho_MMP = moments_MMP.drho, vx1_MMP = moments_MMP.velocityX, vx2_MMP = moments_MMP.velocityY, vx3_MMP = moments_MMP.velocityZ;
    const real drho_PPM = moments_PPM.drho, vx1_PPM = moments_PPM.velocityX, vx2_PPM = moments_PPM.velocityY, vx3_PPM = moments_PPM.velocityZ;
    const real drho_MPM = moments_MPM.drho, vx1_MPM = moments_MPM.velocityX, vx2_MPM = moments_MPM.velocityY, vx3_MPM = moments_MPM.velocityZ;
    const real drho_PMM = moments_PMM.drho, vx1_PMM = moments_PMM.velocityX, vx2_PMM = moments_PMM.velocityY, vx3_PMM = moments_PMM.velocityZ;
    const real drho_MMM = moments_MMM.drho, vx1_MMM = moments_MMM.velocityX, vx2_MMM = moments_MMM.velocityY, vx3_MMM = moments_MMM.velocityZ;

    // second order moments at the source nodes
    const real kxyFromfcNEQ_PPP = moments_PPP.kxyFromfcNEQ, kyzFromfcNEQ_PPP = moments_PPP.kyzFromfcNEQ, kxzFromfcNEQ_PPP = moments_PPP.kxzFromfcNEQ, kxxMyyFromfcNEQ_PPP = moments_PPP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_PPP = moments_PPP.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_MPP = moments_MPP.kxyFromfcNEQ, kyzFromfcNEQ_MPP = moments_MPP.kyzFromfcNEQ, kxzFromfcNEQ_MPP = moments_MPP.kxzFromfcNEQ, kxxMyyFromfcNEQ_MPP = moments_MPP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_MPP = moments_MPP.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_PMP = moments_PMP.kxyFromfcNEQ, kyzFromfcNEQ_PMP = moments_PMP.kyzFromfcNEQ, kxzFromfcNEQ_PMP = moments_PMP.kxzFromfcNEQ, kxxMyyFromfcNEQ_PMP = moments_PMP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_PMP = moments_PMP.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_MMP = moments_MMP.kxyFromfcNEQ, kyzFromfcNEQ_MMP = moments_MMP.kyzFromfcNEQ, kxzFromfcNEQ_MMP = moments_MMP.kxzFromfcNEQ, kxxMyyFromfcNEQ_MMP = moments_MMP.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_MMP = moments_MMP.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_PPM = moments_PPM.kxyFromfcNEQ, kyzFromfcNEQ_PPM = moments_PPM.kyzFromfcNEQ, kxzFromfcNEQ_PPM = moments_PPM.kxzFromfcNEQ, kxxMyyFromfcNEQ_PPM = moments_PPM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_PPM = moments_PPM.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_MPM = moments_MPM.kxyFromfcNEQ, kyzFromfcNEQ_MPM = moments_MPM.kyzFromfcNEQ, kxzFromfcNEQ_MPM = moments_MPM.kxzFromfcNEQ, kxxMyyFromfcNEQ_MPM = moments_MPM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_MPM = moments_MPM.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_PMM = moments_PMM.kxyFromfcNEQ, kyzFromfcNEQ_PMM = moments_PMM.kyzFromfcNEQ, kxzFromfcNEQ_PMM = moments_PMM.kxzFromfcNEQ, kxxMyyFromfcNEQ_PMM = moments_PMM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_PMM = moments_PMM.kxxMzzFromfcNEQ;
    const real kxyFromfcNEQ_MMM = moments_MMM.kxyFromfcNEQ, kyzFromfcNEQ_MMM = moments_MMM.kyzFromfcNEQ, kxzFromfcNEQ_MMM = moments_MMM.kxzFromfcNEQ, kxxMyyFromfcNEQ_MMM = moments_MMM.kxxMyyFromfcNEQ, kxxMzzFromfcNEQ_MMM = moments_MMM.kxxMzzFromfcNEQ;

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

    a_111 = ((vx1_PPP - vx1_MMM) + (vx1_MMP - vx1_PPM)) + ((vx1_MPM - vx1_PMP) + (vx1_PMM - vx1_MPP));
    b_111 = ((vx2_PPP - vx2_MMM) + (vx2_MMP - vx2_PPM)) + ((vx2_MPM - vx2_PMP) + (vx2_PMM - vx2_MPP));
    c_111 = ((vx3_PPP - vx3_MMM) + (vx3_MMP - vx3_PPM)) + ((vx3_MPM - vx3_PMP) + (vx3_PMM - vx3_MPP));

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!- Calculate coefficients for the polynomial interpolation of the pressure
    //! 
    LaplaceRho = 
        ((xoff != c0o1) || (yoff != c0o1) || (zoff != c0o1))
        ? c0o1 : -c3o1 * (a_100 * a_100 + b_010 * b_010 + c_001 * c_001) - c6o1 * (b_100 * a_010 + c_100 * a_001 + c_010 * b_001);
    d_000 = c1o8 * (((drho_PPP + drho_MMM) + (drho_PPM + drho_MMP)) + ((drho_PMM + drho_MPP) + (drho_PMP + drho_MPM)));
    d_100 = c1o4 * (((drho_PPP - drho_MMM) + (drho_PPM - drho_MMP)) + ((drho_PMM - drho_MPP) + (drho_PMP - drho_MPM)));
    d_010 = c1o4 * (((drho_PPP - drho_MMM) + (drho_PPM - drho_MMP)) + ((drho_MPP - drho_PMM) + (drho_MPM - drho_PMP)));
    d_001 = c1o4 * (((drho_PPP - drho_MMM) + (drho_MMP - drho_PPM)) + ((drho_MPP - drho_PMM) + (drho_PMP - drho_MPM)));
    d_110 = c1o2 * (((drho_PPP + drho_MMM) + (drho_PPM + drho_MMP)) - ((drho_PMM + drho_MPP) + (drho_PMP + drho_MPM)));
    d_101 = c1o2 * (((drho_PPP + drho_MMM) - (drho_PPM + drho_MMP)) + ((drho_PMP + drho_MPM) - (drho_PMM + drho_MPP)));
    d_011 = c1o2 * (((drho_PPP + drho_MMM) - (drho_PPM + drho_MMP)) + ((drho_PMM + drho_MPP) - (drho_PMP + drho_MPM)));

    d_111 = (((drho_PPP - drho_MMM) + (drho_MMP - drho_PPM)) + ((drho_PMM - drho_MPP) + (drho_MPM - drho_PMP)));

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

    d_100 = d_100 + yoff * d_110 + zoff * d_101;
    d_010 = d_010 + xoff * d_110 + zoff * d_011;
    d_001 = d_001 + xoff * d_101 + yoff * d_011;
}

}

#endif
