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
//! \file Cumulant27chimStream.cu
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================
/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>
#include "Kernel/Utilities/DistributionHelper.cuh"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
#include "Kernel/ChimeraTransformation.h"

////////////////////////////////////////////////////////////////////////////////
__global__ void LB_Kernel_CumulantK17CompChimRedesigned(
    real omega,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    real* distributions,
    unsigned long numberOfLBnodes,
    int level,
    real* forces,
    real* quadricLimiters,
    real* rho,
    real* veloX,
    real* veloY,
    real* veloZ,
    bool isEvenTimestep,
    const uint *fluidNodeIndices,
    uint numberOfFluidNodes)
{
    //////////////////////////////////////////////////////////////////////////
    //! Cumulant K17 Kernel is based on \ref
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a> and \ref <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
    //!
    //! The cumulant kernel is executed in the following steps
    //!
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get the thread index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned kThread = vf::gpu::getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    // run for all indices in fluidNodeIndices
    if (kThread < numberOfFluidNodes) {
        ////////////////////////////////////////////////////////////////////////////////
        //! - Get the node index from the array containing all indices of fluid nodes
        //!
        const unsigned k = fluidNodeIndices[kThread];

        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist = vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);
        
        ////////////////////////////////////////////////////////////////////////////////
        //! - Set neighbor indices (necessary for indirect addressing)
        //!
        uint k_M00 = neighborX[k];
        uint k_0M0 = neighborY[k];
        uint k_00M = neighborZ[k];
        uint k_MM0 = neighborY[k_M00];
        uint k_M0M = neighborZ[k_M00];
        uint k_0MM = neighborZ[k_0M0];
        uint k_MMM = neighborZ[k_MM0];

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Set local distributions (f's):
        //!
        real f_000 = (dist.f[REST])[k];
        real f_P00 = (dist.f[E])[k];
        real f_M00 = (dist.f[W])[k_M00];
        real f_0P0 = (dist.f[N])[k];
        real f_0M0 = (dist.f[S])[k_0M0];
        real f_00P = (dist.f[T])[k];
        real f_00M = (dist.f[B])[k_00M];
        real f_PP0 = (dist.f[NE])[k];
        real f_MM0 = (dist.f[SW])[k_MM0];
        real f_PM0 = (dist.f[SE])[k_0M0];
        real f_MP0 = (dist.f[NW])[k_M00];
        real f_P0P = (dist.f[TE])[k];
        real f_M0M = (dist.f[BW])[k_M0M];
        real f_P0M = (dist.f[BE])[k_00M];
        real f_M0P = (dist.f[TW])[k_M00];
        real f_0PP = (dist.f[TN])[k];
        real f_0MM = (dist.f[BS])[k_0MM];
        real f_0PM = (dist.f[BN])[k_00M];
        real f_0MP = (dist.f[TS])[k_0M0];
        real f_PPP = (dist.f[TNE])[k];
        real f_MPP = (dist.f[TNW])[k_M00];
        real f_PMP = (dist.f[TSE])[k_0M0];
        real f_MMP = (dist.f[TSW])[k_MM0];
        real f_PPM = (dist.f[BNE])[k_00M];
        real f_MPM = (dist.f[BNW])[k_M0M];
        real f_PMM = (dist.f[BSE])[k_0MM];
        real f_MMM = (dist.f[BSW])[k_MMM];

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Define aliases to use the same variable for the moments (m's):
        //!
        real& m_111 = f_000;
        real& m_211 = f_P00;
        real& m_011 = f_M00;
        real& m_121 = f_0P0;
        real& m_101 = f_0M0;
        real& m_112 = f_00P;
        real& m_110 = f_00M;
        real& m_221 = f_PP0;
        real& m_001 = f_MM0;
        real& m_201 = f_PM0;
        real& m_021 = f_MP0;
        real& m_212 = f_P0P;
        real& m_010 = f_M0M;
        real& m_210 = f_P0M;
        real& m_012 = f_M0P;
        real& m_122 = f_0PP;
        real& m_100 = f_0MM;
        real& m_120 = f_0PM;
        real& m_102 = f_0MP;
        real& m_222 = f_PPP;
        real& m_022 = f_MPP;
        real& m_202 = f_PMP;
        real& m_002 = f_MMP;
        real& m_220 = f_PPM;
        real& m_020 = f_MPM;
        real& m_200 = f_PMM;
        real& m_000 = f_MMM;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
        //!
        real drho = ((((f_PPP + f_MMM) + (f_MPM + f_PMP)) + ((f_MPP + f_PMM) + (f_MMP + f_PPM))) +
                     (((f_0MP + f_0PM) + (f_0MM + f_0PP)) + ((f_M0P + f_P0M) + (f_M0M + f_P0P)) +
                      ((f_MP0 + f_PM0) + (f_MM0 + f_PP0))) +
                     ((f_M00 + f_P00) + (f_0M0 + f_0P0) + (f_00M + f_00P))) +
                    f_000;

        real rho   = c1o1 + drho;
        real OOrho = c1o1 / rho;

        real vvx = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_PMM - f_MPP) + (f_PPM - f_MMP))) +
                    (((f_P0M - f_M0P) + (f_P0P - f_M0M)) + ((f_PM0 - f_MP0) + (f_PP0 - f_MM0))) + (f_P00 - f_M00)) *
                   OOrho;
        real vvy = ((((f_PPP - f_MMM) + (f_MPM - f_PMP)) + ((f_MPP - f_PMM) + (f_PPM - f_MMP))) +
                    (((f_0PM - f_0MP) + (f_0PP - f_0MM)) + ((f_MP0 - f_PM0) + (f_PP0 - f_MM0))) + (f_0P0 - f_0M0)) *
                   OOrho;
        real vvz = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_MPP - f_PMM) + (f_MMP - f_PPM))) +
                    (((f_0MP - f_0PM) + (f_0PP - f_0MM)) + ((f_M0P - f_P0M) + (f_P0P - f_M0M))) + (f_00P - f_00M)) *
                   OOrho;
        ////////////////////////////////////////////////////////////////////////////////////
        //! - Add half of the acceleration (body force) to the velocity as in Eq. (42) \ref
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
        //!
        real factor = c1o1;
        for (size_t i = 1; i <= level; i++) {
            factor *= c2o1;
        }
        real fx = forces[0] / factor;
        real fy = forces[1] / factor;
        real fz = forces[2] / factor;
        vvx += fx * c1o2;
        vvy += fy * c1o2;
        vvz += fz * c1o2;
        ////////////////////////////////////////////////////////////////////////////////////
        // calculate the square of velocities for this lattice node
        real vx2 = vvx * vvx;
        real vy2 = vvy * vvy;
        real vz2 = vvz * vvz;
        ////////////////////////////////////////////////////////////////////////////////////
        //! - Set relaxation limiters for third order cumulants to default value \f$ \lambda=0.001 \f$ according to
        //! section 6 in \ref <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        real wadjust;
        real qudricLimitP = quadricLimiters[0];
        real qudricLimitM = quadricLimiters[1];
        real qudricLimitD = quadricLimiters[2];
        ////////////////////////////////////////////////////////////////////////////////////
        //! - Chimera transform from well conditioned distributions to central moments as defined in Appendix J in \ref
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (6)-(14) in \ref <a
        //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
        //! ]</b></a>
        //!
        ////////////////////////////////////////////////////////////////////////////////////
        // Z - Dir
        forwardInverseChimeraWithK(f_MMM, f_MM0, f_MMP, vvz, vz2, c36o1, c1o36);
        forwardInverseChimeraWithK(f_M0M, f_M00, f_M0P, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(f_MPM, f_MP0, f_MPP, vvz, vz2, c36o1, c1o36);
        forwardInverseChimeraWithK(f_0MM, f_0M0, f_0MP, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(f_00M, f_000, f_00P, vvz, vz2, c9o4, c4o9);
        forwardInverseChimeraWithK(f_0PM, f_0P0, f_0PP, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(f_PMM, f_PM0, f_PMP, vvz, vz2, c36o1, c1o36);
        forwardInverseChimeraWithK(f_P0M, f_P00, f_P0P, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(f_PPM, f_PP0, f_PPP, vvz, vz2, c36o1, c1o36);

        ////////////////////////////////////////////////////////////////////////////////////
        // Y - Dir
        forwardInverseChimeraWithK(f_MMM, f_M0M, f_MPM, vvy, vy2, c6o1, c1o6);
        forwardChimera(f_MM0, f_M00, f_MP0, vvy, vy2);
        forwardInverseChimeraWithK(f_MMP, f_M0P, f_MPP, vvy, vy2, c18o1, c1o18);
        forwardInverseChimeraWithK(f_0MM, f_00M, f_0PM, vvy, vy2, c3o2, c2o3);
        forwardChimera(f_0M0, f_000, f_0P0, vvy, vy2);
        forwardInverseChimeraWithK(f_0MP, f_00P, f_0PP, vvy, vy2, c9o2, c2o9);
        forwardInverseChimeraWithK(f_PMM, f_P0M, f_PPM, vvy, vy2, c6o1, c1o6);
        forwardChimera(f_PM0, f_P00, f_PP0, vvy, vy2);
        forwardInverseChimeraWithK(f_PMP, f_P0P, f_PPP, vvy, vy2, c18o1, c1o18);

        ////////////////////////////////////////////////////////////////////////////////////
        // X - Dir
        forwardInverseChimeraWithK(f_MMM, f_0MM, f_PMM, vvx, vx2, c1o1, c1o1);
        forwardChimera(f_M0M, f_00M, f_P0M, vvx, vx2);
        forwardInverseChimeraWithK(f_MPM, f_0PM, f_PPM, vvx, vx2, c3o1, c1o3);
        forwardChimera(f_MM0, f_0M0, f_PM0, vvx, vx2);
        forwardChimera(f_M00, f_000, f_P00, vvx, vx2);
        forwardChimera(f_MP0, f_0P0, f_PP0, vvx, vx2);
        forwardInverseChimeraWithK(f_MMP, f_0MP, f_PMP, vvx, vx2, c3o1, c1o3);
        forwardChimera(f_M0P, f_00P, f_P0P, vvx, vx2);
        forwardInverseChimeraWithK(f_MPP, f_0PP, f_PPP, vvx, vx2, c3o1, c1o9);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Setting relaxation rates for non-hydrodynamic cumulants (default values). Variable names and equations
        //! according to <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!  => [NAME IN PAPER]=[NAME IN CODE]=[DEFAULT VALUE].
        //!  - Trace of second order cumulants \f$ C_{200}+C_{020}+C_{002} \f$ used to adjust bulk
        //!  viscosity:\f$\omega_2=OxxPyyPzz=1.0 \f$.
        //!  - Third order cumulants \f$ C_{120}+C_{102}, C_{210}+C_{012}, C_{201}+C_{021} \f$: \f$ \omega_3=OxyyPxzz
        //!  \f$ set according to Eq. (111) with simplifications assuming \f$ \omega_2=1.0\f$.
        //!  - Third order cumulants \f$ C_{120}-C_{102}, C_{210}-C_{012}, C_{201}-C_{021} \f$: \f$ \omega_4 = OxyyMxzz
        //!  \f$ set according to Eq. (112) with simplifications assuming \f$ \omega_2 = 1.0\f$.
        //!  - Third order cumulants \f$ C_{111} \f$: \f$ \omega_5 = Oxyz \f$ set according to Eq. (113) with
        //!  simplifications assuming \f$ \omega_2 = 1.0\f$  (modify for different bulk viscosity).
        //!  - Fourth order cumulants \f$ C_{220}, C_{202}, C_{022}, C_{211}, C_{121}, C_{112} \f$: for simplification
        //!  all set to the same default value \f$ \omega_6=\omega_7=\omega_8=O4=1.0 \f$.
        //!  - Fifth order cumulants \f$ C_{221}, C_{212}, C_{122}\f$: \f$\omega_9=O5=1.0\f$.
        //!  - Sixth order cumulant \f$ C_{222}\f$: \f$\omega_{10}=O6=1.0\f$.
        //!
        ////////////////////////////////////////////////////////////
        // 2.
        real OxxPyyPzz = c1o1;
        ////////////////////////////////////////////////////////////
        // 3.
        real OxyyPxzz = c8o1 * (-c2o1 + omega) * (c1o1 + c2o1 * omega) / (-c8o1 - c14o1 * omega + c7o1 * omega * omega);
        real OxyyMxzz =
            c8o1 * (-c2o1 + omega) * (-c7o1 + c4o1 * omega) / (c56o1 - c50o1 * omega + c9o1 * omega * omega);
        real Oxyz = c24o1 * (-c2o1 + omega) * (-c2o1 - c7o1 * omega + c3o1 * omega * omega) /
                    (c48o1 + c152o1 * omega - c130o1 * omega * omega + c29o1 * omega * omega * omega);
        ////////////////////////////////////////////////////////////
        // 4.
        real O4 = c1o1;
        ////////////////////////////////////////////////////////////
        // 5.
        real O5 = c1o1;
        ////////////////////////////////////////////////////////////
        // 6.
        real O6 = c1o1;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - A and B: parameters for fourth order convergence of the diffusion term according to Eq. (115) and (116)
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a> with simplifications assuming \f$ \omega_2 = 1.0 \f$ (modify for
        //! different bulk viscosity).
        //!
        real factorA = (c4o1 + c2o1 * omega - c3o1 * omega * omega) / (c2o1 - c7o1 * omega + c5o1 * omega * omega);
        real factorB = (c4o1 + c28o1 * omega - c14o1 * omega * omega) / (c6o1 - c21o1 * omega + c15o1 * omega * omega);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Compute cumulants from central moments according to Eq. (20)-(23) in
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        ////////////////////////////////////////////////////////////
        // 4.
        real CUMcbb = f_P00 - ((f_PMM + c1o3) * f_M00 + c2o1 * f_00M * f_0M0) * OOrho;
        real CUMbcb = f_0P0 - ((f_MPM + c1o3) * f_0M0 + c2o1 * f_00M * f_M00) * OOrho;
        real CUMbbc = f_00P - ((f_MMP + c1o3) * f_00M + c2o1 * f_0M0 * f_M00) * OOrho;

        real CUMcca =
            f_PPM - (((f_PMM * f_MPM + c2o1 * f_00M * f_00M) + c1o3 * (f_PMM + f_MPM)) * OOrho - c1o9 * (drho * OOrho));
        real CUMcac =
            f_PMP - (((f_PMM * f_MMP + c2o1 * f_0M0 * f_0M0) + c1o3 * (f_PMM + f_MMP)) * OOrho - c1o9 * (drho * OOrho));
        real CUMacc =
            f_MPP - (((f_MMP * f_MPM + c2o1 * f_M00 * f_M00) + c1o3 * (f_MMP + f_MPM)) * OOrho - c1o9 * (drho * OOrho));
        ////////////////////////////////////////////////////////////
        // 5.
        real CUMbcc =
            f_0PP - ((f_MMP * f_0PM + f_MPM * f_0MP + c4o1 * f_M00 * f_000 + c2o1 * (f_0M0 * f_MP0 + f_00M * f_M0P)) +
                     c1o3 * (f_0PM + f_0MP)) *
                        OOrho;
        real CUMcbc =
            f_P0P - ((f_MMP * f_P0M + f_PMM * f_M0P + c4o1 * f_0M0 * f_000 + c2o1 * (f_M00 * f_PM0 + f_00M * f_0MP)) +
                     c1o3 * (f_P0M + f_M0P)) *
                        OOrho;
        real CUMccb =
            f_PP0 - ((f_PMM * f_MP0 + f_MPM * f_PM0 + c4o1 * f_00M * f_000 + c2o1 * (f_0M0 * f_0PM + f_M00 * f_P0M)) +
                     c1o3 * (f_MP0 + f_PM0)) *
                        OOrho;
        ////////////////////////////////////////////////////////////
        // 6.
        real CUMccc = f_PPP + ((-c4o1 * f_000 * f_000 - (f_PMM * f_MPP + f_MPM * f_PMP + f_MMP * f_PPM) -
                                c4o1 * (f_M00 * f_P00 + f_0M0 * f_0P0 + f_00M * f_00P) -
                                c2o1 * (f_0PM * f_0MP + f_P0M * f_M0P + f_PM0 * f_MP0)) *
                                   OOrho +
                               (c4o1 * (f_0M0 * f_0M0 * f_MPM + f_M00 * f_M00 * f_PMM + f_00M * f_00M * f_MMP) +
                                c2o1 * (f_PMM * f_MPM * f_MMP) + c16o1 * f_00M * f_0M0 * f_M00) *
                                   OOrho * OOrho -
                               c1o3 * (f_MPP + f_PMP + f_PPM) * OOrho - c1o9 * (f_PMM + f_MPM + f_MMP) * OOrho +
                               (c2o1 * (f_0M0 * f_0M0 + f_M00 * f_M00 + f_00M * f_00M) +
                                (f_MMP * f_MPM + f_MMP * f_PMM + f_MPM * f_PMM) + c1o3 * (f_MMP + f_MPM + f_PMM)) *
                                   OOrho * OOrho * c2o3 +
                               c1o27 * ((drho * drho - drho) * OOrho * OOrho));

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Compute linear combinations of second and third order cumulants
        //!
        ////////////////////////////////////////////////////////////
        // 2.
        real mxxPyyPzz = f_PMM + f_MPM + f_MMP;
        real mxxMyy    = f_PMM - f_MPM;
        real mxxMzz    = f_PMM - f_MMP;
        ////////////////////////////////////////////////////////////
        // 3.
        real mxxyPyzz = f_P0M + f_M0P;
        real mxxyMyzz = f_P0M - f_M0P;

        real mxxzPyyz = f_PM0 + f_MP0;
        real mxxzMyyz = f_PM0 - f_MP0;

        real mxyyPxzz = f_0PM + f_0MP;
        real mxyyMxzz = f_0PM - f_0MP;

        ////////////////////////////////////////////////////////////////////////////////////
        // incl. correction
        ////////////////////////////////////////////////////////////
        //! - Compute velocity  gradients from second order cumulants according to Eq. (27)-(32)
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a> Further explanations of the correction in viscosity in Appendix H of
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> Note that the division by rho is omitted here as we need rho times
        //! the gradients later.
        //!
        real Dxy  = -c3o1 * omega * f_00M;
        real Dxz  = -c3o1 * omega * f_0M0;
        real Dyz  = -c3o1 * omega * f_M00;
        real dxux = c1o2 * (-omega) * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (f_MMM - mxxPyyPzz);
        real dyuy = dxux + omega * c3o2 * mxxMyy;
        real dzuz = dxux + omega * c3o2 * mxxMzz;
        ////////////////////////////////////////////////////////////
        //! - Relaxation of second order cumulants with correction terms according to Eq. (33)-(35) in
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        mxxPyyPzz +=
            OxxPyyPzz * (f_MMM - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
        mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
        mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

        ////////////////////////////////////////////////////////////////////////////////////
        ////no correction
        // mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);
        // mxxMyy += -(-omega) * (-mxxMyy);
        // mxxMzz += -(-omega) * (-mxxMzz);
        //////////////////////////////////////////////////////////////////////////
        f_M00 += omega * (-f_M00);
        f_0M0 += omega * (-f_0M0);
        f_00M += omega * (-f_00M);

        ////////////////////////////////////////////////////////////////////////////////////
        // relax
        //////////////////////////////////////////////////////////////////////////
        // incl. limiter
        //! - Relaxation of third order cumulants including limiter according to Eq. (116)-(123)
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        wadjust = Oxyz + (c1o1 - Oxyz) * abs(f_000) / (abs(f_000) + qudricLimitD);
        f_000 += wadjust * (-f_000);
        wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxxyPyzz) / (abs(mxxyPyzz) + qudricLimitP);
        mxxyPyzz += wadjust * (-mxxyPyzz);
        wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxxyMyzz) / (abs(mxxyMyzz) + qudricLimitM);
        mxxyMyzz += wadjust * (-mxxyMyzz);
        wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxxzPyyz) / (abs(mxxzPyyz) + qudricLimitP);
        mxxzPyyz += wadjust * (-mxxzPyyz);
        wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxxzMyyz) / (abs(mxxzMyyz) + qudricLimitM);
        mxxzMyyz += wadjust * (-mxxzMyyz);
        wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxyyPxzz) / (abs(mxyyPxzz) + qudricLimitP);
        mxyyPxzz += wadjust * (-mxyyPxzz);
        wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxyyMxzz) / (abs(mxyyMxzz) + qudricLimitM);
        mxyyMxzz += wadjust * (-mxyyMxzz);
        //////////////////////////////////////////////////////////////////////////
        // no limiter
        // mfbbb += OxyyMxzz * (-mfbbb);
        // mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
        // mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
        // mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
        // mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
        // mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
        // mxyyMxzz += OxyyMxzz * (-mxyyMxzz);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Compute inverse linear combinations of second and third order cumulants
        //!
        f_PMM = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
        f_MPM = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz);
        f_MMP = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

        f_P0M = (mxxyMyzz + mxxyPyzz) * c1o2;
        f_M0P = (-mxxyMyzz + mxxyPyzz) * c1o2;
        f_PM0 = (mxxzMyyz + mxxzPyyz) * c1o2;
        f_MP0 = (-mxxzMyyz + mxxzPyyz) * c1o2;
        f_0PM = (mxyyMxzz + mxyyPxzz) * c1o2;
        f_0MP = (-mxyyMxzz + mxyyPxzz) * c1o2;
        //////////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////////
        // 4.
        // no limiter
        //! - Relax fourth order cumulants to modified equilibrium for fourth order convergence of diffusion according
        //! to Eq. (43)-(48) <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        CUMacc = -O4 * (c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * factorA + (c1o1 - O4) * (CUMacc);
        CUMcac = -O4 * (c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * factorA + (c1o1 - O4) * (CUMcac);
        CUMcca = -O4 * (c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * factorA + (c1o1 - O4) * (CUMcca);
        CUMbbc = -O4 * (c1o1 / omega - c1o2) * Dxy * c1o3 * factorB + (c1o1 - O4) * (CUMbbc);
        CUMbcb = -O4 * (c1o1 / omega - c1o2) * Dxz * c1o3 * factorB + (c1o1 - O4) * (CUMbcb);
        CUMcbb = -O4 * (c1o1 / omega - c1o2) * Dyz * c1o3 * factorB + (c1o1 - O4) * (CUMcbb);

        //////////////////////////////////////////////////////////////////////////
        // 5.
        CUMbcc += O5 * (-CUMbcc);
        CUMcbc += O5 * (-CUMcbc);
        CUMccb += O5 * (-CUMccb);

        //////////////////////////////////////////////////////////////////////////
        // 6.
        CUMccc += O6 * (-CUMccc);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Compute central moments from post collision cumulants according to Eq. (53)-(56) in
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!

        //////////////////////////////////////////////////////////////////////////
        // 4.
        f_P00 = CUMcbb + c1o3 * ((c3o1 * f_PMM + c1o1) * f_M00 + c6o1 * f_00M * f_0M0) * OOrho;
        f_0P0 = CUMbcb + c1o3 * ((c3o1 * f_MPM + c1o1) * f_0M0 + c6o1 * f_00M * f_M00) * OOrho;
        f_00P = CUMbbc + c1o3 * ((c3o1 * f_MMP + c1o1) * f_00M + c6o1 * f_0M0 * f_M00) * OOrho;

        f_PPM =
            CUMcca +
            (((f_PMM * f_MPM + c2o1 * f_00M * f_00M) * c9o1 + c3o1 * (f_PMM + f_MPM)) * OOrho - (drho * OOrho)) * c1o9;
        f_PMP =
            CUMcac +
            (((f_PMM * f_MMP + c2o1 * f_0M0 * f_0M0) * c9o1 + c3o1 * (f_PMM + f_MMP)) * OOrho - (drho * OOrho)) * c1o9;
        f_MPP =
            CUMacc +
            (((f_MMP * f_MPM + c2o1 * f_M00 * f_M00) * c9o1 + c3o1 * (f_MMP + f_MPM)) * OOrho - (drho * OOrho)) * c1o9;

        //////////////////////////////////////////////////////////////////////////
        // 5.
        f_0PP = CUMbcc + c1o3 *
                             (c3o1 * (f_MMP * f_0PM + f_MPM * f_0MP + c4o1 * f_M00 * f_000 +
                                      c2o1 * (f_0M0 * f_MP0 + f_00M * f_M0P)) +
                              (f_0PM + f_0MP)) *
                             OOrho;
        f_P0P = CUMcbc + c1o3 *
                             (c3o1 * (f_MMP * f_P0M + f_PMM * f_M0P + c4o1 * f_0M0 * f_000 +
                                      c2o1 * (f_M00 * f_PM0 + f_00M * f_0MP)) +
                              (f_P0M + f_M0P)) *
                             OOrho;
        f_PP0 = CUMccb + c1o3 *
                             (c3o1 * (f_PMM * f_MP0 + f_MPM * f_PM0 + c4o1 * f_00M * f_000 +
                                      c2o1 * (f_0M0 * f_0PM + f_M00 * f_P0M)) +
                              (f_MP0 + f_PM0)) *
                             OOrho;

        //////////////////////////////////////////////////////////////////////////
        // 6.
        f_PPP = CUMccc - ((-c4o1 * f_000 * f_000 - (f_PMM * f_MPP + f_MPM * f_PMP + f_MMP * f_PPM) -
                           c4o1 * (f_M00 * f_P00 + f_0M0 * f_0P0 + f_00M * f_00P) -
                           c2o1 * (f_0PM * f_0MP + f_P0M * f_M0P + f_PM0 * f_MP0)) *
                              OOrho +
                          (c4o1 * (f_0M0 * f_0M0 * f_MPM + f_M00 * f_M00 * f_PMM + f_00M * f_00M * f_MMP) +
                           c2o1 * (f_PMM * f_MPM * f_MMP) + c16o1 * f_00M * f_0M0 * f_M00) *
                              OOrho * OOrho -
                          c1o3 * (f_MPP + f_PMP + f_PPM) * OOrho - c1o9 * (f_PMM + f_MPM + f_MMP) * OOrho +
                          (c2o1 * (f_0M0 * f_0M0 + f_M00 * f_M00 + f_00M * f_00M) +
                           (f_MMP * f_MPM + f_MMP * f_PMM + f_MPM * f_PMM) + c1o3 * (f_MMP + f_MPM + f_PMM)) *
                              OOrho * OOrho * c2o3 +
                          c1o27 * ((drho * drho - drho) * OOrho * OOrho));

        ////////////////////////////////////////////////////////////////////////////////////
        //! -  Add acceleration (body force) to first order cumulants according to Eq. (85)-(87) in
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
        //!
        f_0MM = -f_0MM;
        f_M0M = -f_M0M;
        f_MM0 = -f_MM0;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
        //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
        //! ]</b></a>
        //!
        ////////////////////////////////////////////////////////////////////////////////////
        // X - Dir
        backwardInverseChimeraWithK(f_MMM, f_0MM, f_PMM, vvx, vx2, c1o1, c1o1);
        backwardChimera(f_M0M, f_00M, f_P0M, vvx, vx2);
        backwardInverseChimeraWithK(f_MPM, f_0PM, f_PPM, vvx, vx2, c3o1, c1o3);
        backwardChimera(f_MM0, f_0M0, f_PM0, vvx, vx2);
        backwardChimera(f_M00, f_000, f_P00, vvx, vx2);
        backwardChimera(f_MP0, f_0P0, f_PP0, vvx, vx2);
        backwardInverseChimeraWithK(f_MMP, f_0MP, f_PMP, vvx, vx2, c3o1, c1o3);
        backwardChimera(f_M0P, f_00P, f_P0P, vvx, vx2);
        backwardInverseChimeraWithK(f_MPP, f_0PP, f_PPP, vvx, vx2, c9o1, c1o9);

        ////////////////////////////////////////////////////////////////////////////////////
        // Y - Dir
        backwardInverseChimeraWithK(f_MMM, f_M0M, f_MPM, vvy, vy2, c6o1, c1o6);
        backwardChimera(f_MM0, f_M00, f_MP0, vvy, vy2);
        backwardInverseChimeraWithK(f_MMP, f_M0P, f_MPP, vvy, vy2, c18o1, c1o18);
        backwardInverseChimeraWithK(f_0MM, f_00M, f_0PM, vvy, vy2, c3o2, c2o3);
        backwardChimera(f_0M0, f_000, f_0P0, vvy, vy2);
        backwardInverseChimeraWithK(f_0MP, f_00P, f_0PP, vvy, vy2, c9o2, c2o9);
        backwardInverseChimeraWithK(f_PMM, f_P0M, f_PPM, vvy, vy2, c6o1, c1o6);
        backwardChimera(f_PM0, f_P00, f_PP0, vvy, vy2);
        backwardInverseChimeraWithK(f_PMP, f_P0P, f_PPP, vvy, vy2, c18o1, c1o18);

        ////////////////////////////////////////////////////////////////////////////////////
        // Z - Dir
        backwardInverseChimeraWithK(f_MMM, f_MM0, f_MMP, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(f_M0M, f_M00, f_M0P, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(f_MPM, f_MP0, f_MPP, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(f_0MM, f_0M0, f_0MP, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(f_00M, f_000, f_00P, vvz, vz2, c9o4, c4o9);
        backwardInverseChimeraWithK(f_0PM, f_0P0, f_0PP, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(f_PMM, f_PM0, f_PMP, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(f_P0M, f_P00, f_P0P, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(f_PPM, f_PP0, f_PPP, vvz, vz2, c36o1, c1o36);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Write distributions: style of reading and writing the distributions from/to
        //! stored arrays dependent on timestep is based on the esoteric twist algorithm
        //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        (dist.f[E])[k]      = f_M00;
        (dist.f[W])[k_M00]     = f_P00;
        (dist.f[N])[k]      = f_0M0;
        (dist.f[S])[k_0M0]     = f_0P0;
        (dist.f[T])[k]      = f_00M;
        (dist.f[B])[k_00M]     = f_00P;
        (dist.f[NE])[k]     = f_MM0;
        (dist.f[SW])[k_MM0]   = f_PP0;
        (dist.f[SE])[k_0M0]    = f_MP0;
        (dist.f[NW])[k_M00]    = f_PM0;
        (dist.f[TE])[k]     = f_M0M;
        (dist.f[BW])[k_M0M]   = f_P0P;
        (dist.f[BE])[k_00M]    = f_M0P;
        (dist.f[TW])[k_M00]    = f_P0M;
        (dist.f[TN])[k]     = f_0MM;
        (dist.f[BS])[k_0MM]   = f_0PP;
        (dist.f[BN])[k_00M]    = f_0MP;
        (dist.f[TS])[k_0M0]    = f_0PM;
        (dist.f[REST])[k]   = f_000;
        (dist.f[TNE])[k]    = f_MMM;
        (dist.f[TSE])[k_0M0]   = f_MPM;
        (dist.f[BNE])[k_00M]   = f_MMP;
        (dist.f[BSE])[k_0MM]  = f_MPP;
        (dist.f[TNW])[k_M00]   = f_PMM;
        (dist.f[TSW])[k_MM0]  = f_PPM;
        (dist.f[BNW])[k_M0M]  = f_PMP;
        (dist.f[BSW])[k_MMM] = f_PPP;
    }
}