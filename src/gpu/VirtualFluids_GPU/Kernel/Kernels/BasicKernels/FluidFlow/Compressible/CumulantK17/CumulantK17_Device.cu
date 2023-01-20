
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
//! \file CumlantK17_Device.cu
//! \author Anna Wellmann, Martin Schönherr, Henry Korb, Henrik Asmuth
//! \date 05/12/2022
//! \brief Kernel for CumulantK17 including different turbulence models and options for local body forces and writing macroscopic variables
//!
//! CumulantK17 kernel using chimera transformations and quartic limiters as present in Geier et al. (2017). Additional options are three different
//! eddy-viscosity turbulence models (Smagorinsky, AMD, QR) that can be set via the template parameter turbulenceModel (with default
//! TurbulenceModel::None).
//! The kernel is executed separately for each subset of fluid node indices with a different tag CollisionTemplate. For each subset, only the locally
//! required options are switched on ( \param writeMacroscopicVariables and/or \param applyBodyForce) in order to minimize memory accesses. The default
//! refers to the plain cumlant kernel (CollisionTemplate::Default).
//! Nodes are added to subsets (taggedFluidNodes) in Simulation::init using a corresponding tag with different values of CollisionTemplate. These subsets
//! are provided by the utilized PostCollisionInteractiors depending on they specifc requirements (e.g. writeMacroscopicVariables for probes).

//=======================================================================================
/* Device code */
#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>
#include "Kernel/Utilities/DistributionHelper.cuh"

#include "GPU/TurbulentViscosityInlines.cuh"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
#include "Kernel/Utilities/ChimeraTransformation.h"


////////////////////////////////////////////////////////////////////////////////
template<TurbulenceModel turbulenceModel, bool writeMacroscopicVariables, bool applyBodyForce>
__global__ void LB_Kernel_CumulantK17(
    real omega_in,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    real* distributions,
    real* rho,
    real* vx,
    real* vy,
    real* vz,
    real* turbulentViscosity,
    real SGSconstant,
    unsigned long long numberOfLBnodes,
    int level,
    real* forces,
    real* bodyForceX,
    real* bodyForceY,
    real* bodyForceZ,
    real* quadricLimiters,
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
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned kThread = vf::gpu::getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    // run for all indices in size_Mat and fluid nodes
    if (kThread >= numberOfFluidNodes)
        return;
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get the node index from the array containing all indices of fluid nodes
    //!
    const unsigned k_000 = fluidNodeIndices[kThread];

    //////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
    //! timestep is based on the esoteric twist algorithm \ref <a
    //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    Distributions27 dist = vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set neighbor indices (necessary for indirect addressing)
    uint k_M00 = neighborX[k_000];
    uint k_0M0 = neighborY[k_000];
    uint k_00M = neighborZ[k_000];
    uint k_MM0 = neighborY[k_M00];
    uint k_M0M = neighborZ[k_M00];
    uint k_0MM = neighborZ[k_0M0];
    uint k_MMM = neighborZ[k_MM0];
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set local distributions
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

    //////////////////////////////////////////////////////(unsigned long)//////////////////////////////
    //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    real drho = ((((f_PPP + f_MMM) + (f_MPM + f_PMP)) + ((f_MPP + f_PMM) + (f_MMP + f_PPM))) +
                (((f_0MP + f_0PM) + (f_0MM + f_0PP)) + ((f_M0P + f_P0M) + (f_M0M + f_P0P)) +
                ((f_MP0 + f_PM0) + (f_MM0 + f_PP0))) +
                ((f_M00 + f_P00) + (f_0M0 + f_0P0) + (f_00M + f_00P))) +
                    f_000;

    real oneOverRho = c1o1 / (c1o1 + drho);

    real vvx = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_PMM - f_MPP) + (f_PPM - f_MMP))) +
                (((f_P0M - f_M0P) + (f_P0P - f_M0M)) + ((f_PM0 - f_MP0) + (f_PP0 - f_MM0))) + (f_P00 - f_M00)) *
            oneOverRho;
    real vvy = ((((f_PPP - f_MMM) + (f_MPM - f_PMP)) + ((f_MPP - f_PMM) + (f_PPM - f_MMP))) +
                (((f_0PM - f_0MP) + (f_0PP - f_0MM)) + ((f_MP0 - f_PM0) + (f_PP0 - f_MM0))) + (f_0P0 - f_0M0)) *
            oneOverRho;
    real vvz = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_MPP - f_PMM) + (f_MMP - f_PPM))) +
                (((f_0MP - f_0PM) + (f_0PP - f_0MM)) + ((f_M0P - f_P0M) + (f_P0P - f_M0M))) + (f_00P - f_00M)) *
            oneOverRho;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Add half of the acceleration (body force) to the velocity as in Eq. (42) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    real factor = c1o1;
    for (size_t i = 1; i <= level; i++) {
        factor *= c2o1;
    }

    real fx = forces[0];
    real fy = forces[1];
    real fz = forces[2];

    if( applyBodyForce ){
        fx += bodyForceX[k_000];
        fy += bodyForceY[k_000];
        fz += bodyForceZ[k_000];

        // real vx = vvx;
        // real vy = vvy;
        // real vz = vvz;
        real acc_x = fx * c1o2 / factor;
        real acc_y = fy * c1o2 / factor;
        real acc_z = fz * c1o2 / factor;

        vvx += acc_x;
        vvy += acc_y;
        vvz += acc_z;

        // Reset body force. To be used when not using round-off correction.
        bodyForceX[k_000] = 0.0f;
        bodyForceY[k_000] = 0.0f;
        bodyForceZ[k_000] = 0.0f;

        ////////////////////////////////////////////////////////////////////////////////////
        //!> Round-off correction
        //!
        //!> Similar to Kahan summation algorithm (https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
        //!> Essentially computes the round-off error of the applied force and adds it in the next time step as a compensation.
        //!> Seems to be necesseary at very high Re boundary layers, where the forcing and velocity can
        //!> differ by several orders of magnitude.
        //!> \note 16/05/2022: Testing, still ongoing!
        //!
        // bodyForceX[k_000] = (acc_x-(vvx-vx))*factor*c2o1;
        // bodyForceY[k_000] = (acc_y-(vvy-vy))*factor*c2o1;
        // bodyForceZ[k_000] = (acc_z-(vvz-vz))*factor*c2o1;
    }
    else{
        vvx += fx * c1o2 / factor;
        vvy += fy * c1o2 / factor;
        vvz += fz * c1o2 / factor;
    }


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
    real quadricLimitP = quadricLimiters[0];
    real quadricLimitM = quadricLimiters[1];
    real quadricLimitD = quadricLimiters[2];
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
    forwardInverseChimeraWithK(f_M0M, f_M00, f_M0P, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(f_MPM, f_MP0, f_MPP, vvz, vz2, c36o1, c1o36);
    forwardInverseChimeraWithK(f_0MM, f_0M0, f_0MP, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(f_00M, f_000, f_00P, vvz, vz2, c9o4,  c4o9);
    forwardInverseChimeraWithK(f_0PM, f_0P0, f_0PP, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(f_PMM, f_PM0, f_PMP, vvz, vz2, c36o1, c1o36);
    forwardInverseChimeraWithK(f_P0M, f_P00, f_P0P, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(f_PPM, f_PP0, f_PPP, vvz, vz2, c36o1, c1o36);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    forwardInverseChimeraWithK(f_MMM, f_M0M, f_MPM, vvy, vy2, c6o1,  c1o6);
    forwardChimera(            f_MM0, f_M00, f_MP0, vvy, vy2);
    forwardInverseChimeraWithK(f_MMP, f_M0P, f_MPP, vvy, vy2, c18o1, c1o18);
    forwardInverseChimeraWithK(f_0MM, f_00M, f_0PM, vvy, vy2, c3o2,  c2o3);
    forwardChimera(            f_0M0, f_000, f_0P0, vvy, vy2);
    forwardInverseChimeraWithK(f_0MP, f_00P, f_0PP, vvy, vy2, c9o2,  c2o9);
    forwardInverseChimeraWithK(f_PMM, f_P0M, f_PPM, vvy, vy2, c6o1,  c1o6);
    forwardChimera(            f_PM0, f_P00, f_PP0, vvy, vy2);
    forwardInverseChimeraWithK(f_PMP, f_P0P, f_PPP, vvy, vy2, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    forwardInverseChimeraWithK(f_MMM, f_0MM, f_PMM, vvx, vx2, c1o1, c1o1);
    forwardChimera(            f_M0M, f_00M, f_P0M, vvx, vx2);
    forwardInverseChimeraWithK(f_MPM, f_0PM, f_PPM, vvx, vx2, c3o1, c1o3);
    forwardChimera(            f_MM0, f_0M0, f_PM0, vvx, vx2);
    forwardChimera(            f_M00, f_000, f_P00, vvx, vx2);
    forwardChimera(            f_MP0, f_0P0, f_PP0, vvx, vx2);
    forwardInverseChimeraWithK(f_MMP, f_0MP, f_PMP, vvx, vx2, c3o1, c1o3);
    forwardChimera(            f_M0P, f_00P, f_P0P, vvx, vx2);
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
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Calculate modified omega with turbulent viscosity
    //!
    real omega = omega_in;
    if(turbulenceModel != TurbulenceModel::None){ omega /= (c1o1 + c3o1*omega_in*turbulentViscosity[k_000]); }
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
    //! - A and DIR_00M: parameters for fourth order convergence of the diffusion term according to Eq. (115) and (116)
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
    real c_211 = m_211 - ((m_200 + c1o3) * m_011 + c2o1 * m_110 * m_101) * oneOverRho;
    real c_121 = m_121 - ((m_020 + c1o3) * m_101 + c2o1 * m_110 * m_011) * oneOverRho;
    real c_112 = m_112 - ((m_002 + c1o3) * m_110 + c2o1 * m_101 * m_011) * oneOverRho;

    real c_220 = m_220 - (((m_200 * m_020 + c2o1 * m_110 * m_110) + c1o3 * (m_200 + m_020)) * oneOverRho - c1o9 * (drho * oneOverRho));
    real c_202 = m_202 - (((m_200 * m_002 + c2o1 * m_101 * m_101) + c1o3 * (m_200 + m_002)) * oneOverRho - c1o9 * (drho * oneOverRho));
    real c_022 = m_022 - (((m_002 * m_020 + c2o1 * m_011 * m_011) + c1o3 * (m_002 + m_020)) * oneOverRho - c1o9 * (drho * oneOverRho));
    ////////////////////////////////////////////////////////////
    // 5.
    real c_122 =
        m_122 - ((m_002 * m_120 + m_020 * m_102 + c4o1 * m_011 * m_111 + c2o1 * (m_101 * m_021 + m_110 * m_012)) +
                c1o3 * (m_120 + m_102)) *
                oneOverRho;
    real c_212 =
        m_212 - ((m_002 * m_210 + m_200 * m_012 + c4o1 * m_101 * m_111 + c2o1 * (m_011 * m_201 + m_110 * m_102)) +
                c1o3 * (m_210 + m_012)) *
                oneOverRho;
    real c_221 =
        m_221 - ((m_200 * m_021 + m_020 * m_201 + c4o1 * m_110 * m_111 + c2o1 * (m_101 * m_120 + m_011 * m_210)) +
                c1o3 * (m_021 + m_201)) *
                oneOverRho;
    ////////////////////////////////////////////////////////////
    // 6.
    real c_222 = m_222 + ((-c4o1 * m_111 * m_111 - (m_200 * m_022 + m_020 * m_202 + m_002 * m_220) -
                            c4o1 * (m_011 * m_211 + m_101 * m_121 + m_110 * m_112) -
                            c2o1 * (m_120 * m_102 + m_210 * m_012 + m_201 * m_021)) *
                            oneOverRho +
                        (c4o1 * (m_101 * m_101 * m_020 + m_011 * m_011 * m_200 + m_110 * m_110 * m_002) +
                            c2o1 * (m_200 * m_020 * m_002) + c16o1 * m_110 * m_101 * m_011) *
                            oneOverRho * oneOverRho -
                            c1o3 * (m_022 + m_202 + m_220) * oneOverRho - c1o9 * (m_200 + m_020 + m_002) * oneOverRho +
                        (c2o1 * (m_101 * m_101 + m_011 * m_011 + m_110 * m_110) +
                            (m_002 * m_020 + m_002 * m_200 + m_020 * m_200) + c1o3 * (m_002 + m_020 + m_200)) *
                            oneOverRho * oneOverRho * c2o3 +
                            c1o27 * ((drho * drho - drho) * oneOverRho * oneOverRho));

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute linear combinations of second and third order cumulants
    //!
    ////////////////////////////////////////////////////////////
    // 2.
    real mxxPyyPzz = m_200 + m_020 + m_002;
    real mxxMyy    = m_200 - m_020;
    real mxxMzz    = m_200 - m_002;
    ////////////////////////////////////////////////////////////
    // 3.
    real mxxyPyzz = m_210 + m_012;
    real mxxyMyzz = m_210 - m_012;

    real mxxzPyyz = m_201 + m_021;
    real mxxzMyyz = m_201 - m_021;

    real mxyyPxzz = m_120 + m_102;
    real mxyyMxzz = m_120 - m_102;

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
    real Dxy  = -c3o1 * omega * m_110;
    real Dxz  = -c3o1 * omega * m_101;
    real Dyz  = -c3o1 * omega * m_011;
    real dxux = c1o2 * (-omega) * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (m_000 - mxxPyyPzz);
    real dyuy = dxux + omega * c3o2 * mxxMyy;
    real dzuz = dxux + omega * c3o2 * mxxMzz;

    ////////////////////////////////////////////////////////////////////////////////////
    switch (turbulenceModel)
    {
    case TurbulenceModel::None:
    case TurbulenceModel::AMD:  //AMD is computed in separate kernel
        break;
    case TurbulenceModel::Smagorinsky:
        turbulentViscosity[k_000] = calcTurbulentViscositySmagorinsky(SGSconstant, dxux, dyuy, dzuz, Dxy, Dxz , Dyz);
        break;
    case TurbulenceModel::QR:
        turbulentViscosity[k_000] = calcTurbulentViscosityQR(SGSconstant, dxux, dyuy, dzuz, Dxy, Dxz , Dyz);
        break;
    default:
        break;
    }
    ////////////////////////////////////////////////////////////
    //! - Relaxation of second order cumulants with correction terms according to Eq. (33)-(35) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!
    mxxPyyPzz += OxxPyyPzz * (m_000 - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
    mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
    mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

    ////////////////////////////////////////////////////////////////////////////////////
    ////no correction
    // mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);
    // mxxMyy += -(-omega) * (-mxxMyy);
    // mxxMzz += -(-omega) * (-mxxMzz);
    //////////////////////////////////////////////////////////////////////////
    m_011 += omega * (-m_011);
    m_101 += omega * (-m_101);
    m_110 += omega * (-m_110);

    ////////////////////////////////////////////////////////////////////////////////////
    // relax
    //////////////////////////////////////////////////////////////////////////
    // incl. limiter
    //! - Relaxation of third order cumulants including limiter according to Eq. (116)-(123)
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!
    real wadjust = Oxyz + (c1o1 - Oxyz) * abs(m_111) / (abs(m_111) + quadricLimitD);
    m_111 += wadjust * (-m_111);
    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxxyPyzz) / (abs(mxxyPyzz) + quadricLimitP);
    mxxyPyzz += wadjust * (-mxxyPyzz);
    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxxyMyzz) / (abs(mxxyMyzz) + quadricLimitM);
    mxxyMyzz += wadjust * (-mxxyMyzz);
    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxxzPyyz) / (abs(mxxzPyyz) + quadricLimitP);
    mxxzPyyz += wadjust * (-mxxzPyyz);
    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxxzMyyz) / (abs(mxxzMyyz) + quadricLimitM);
    mxxzMyyz += wadjust * (-mxxzMyyz);
    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxyyPxzz) / (abs(mxyyPxzz) + quadricLimitP);
    mxyyPxzz += wadjust * (-mxyyPxzz);
    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxyyMxzz) / (abs(mxyyMxzz) + quadricLimitM);
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
    m_200 = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
    m_020 = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz);
    m_002 = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

    m_210 = ( mxxyMyzz + mxxyPyzz) * c1o2;
    m_012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
    m_201 = ( mxxzMyyz + mxxzPyyz) * c1o2;
    m_021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
    m_120 = ( mxyyMxzz + mxyyPxzz) * c1o2;
    m_102 = (-mxyyMxzz + mxyyPxzz) * c1o2;
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // 4.
    // no limiter
    //! - Relax fourth order cumulants to modified equilibrium for fourth order convergence of diffusion according
    //! to Eq. (43)-(48) <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!
    c_022 = -O4 * (c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * factorA + (c1o1 - O4) * (c_022);
    c_202 = -O4 * (c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * factorA + (c1o1 - O4) * (c_202);
    c_220 = -O4 * (c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * factorA + (c1o1 - O4) * (c_220);
    c_112 = -O4 * (c1o1 / omega - c1o2) * Dxy           * c1o3 * factorB + (c1o1 - O4) * (c_112);
    c_121 = -O4 * (c1o1 / omega - c1o2) * Dxz           * c1o3 * factorB + (c1o1 - O4) * (c_121);
    c_211 = -O4 * (c1o1 / omega - c1o2) * Dyz           * c1o3 * factorB + (c1o1 - O4) * (c_211);


    //////////////////////////////////////////////////////////////////////////
    // 5.
    c_122 += O5 * (-c_122);
    c_212 += O5 * (-c_212);
    c_221 += O5 * (-c_221);

    //////////////////////////////////////////////////////////////////////////
    // 6.
    c_222 += O6 * (-c_222);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute central moments from post collision cumulants according to Eq. (53)-(56) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!

    //////////////////////////////////////////////////////////////////////////
    // 4.
    m_211 = c_211 + c1o3 * ((c3o1 * m_200 + c1o1) * m_011 + c6o1 * m_110 * m_101) * oneOverRho;
    m_121 = c_121 + c1o3 * ((c3o1 * m_020 + c1o1) * m_101 + c6o1 * m_110 * m_011) * oneOverRho;
    m_112 = c_112 + c1o3 * ((c3o1 * m_002 + c1o1) * m_110 + c6o1 * m_101 * m_011) * oneOverRho;

    m_220 =
        c_220 + (((m_200 * m_020 + c2o1 * m_110 * m_110) * c9o1 + c3o1 * (m_200 + m_020)) * oneOverRho - (drho * oneOverRho)) * c1o9;
    m_202 =
        c_202 + (((m_200 * m_002 + c2o1 * m_101 * m_101) * c9o1 + c3o1 * (m_200 + m_002)) * oneOverRho - (drho * oneOverRho)) * c1o9;
    m_022 =
        c_022 + (((m_002 * m_020 + c2o1 * m_011 * m_011) * c9o1 + c3o1 * (m_002 + m_020)) * oneOverRho - (drho * oneOverRho)) * c1o9;

    //////////////////////////////////////////////////////////////////////////
    // 5.
    m_122 = c_122 + c1o3 *
            (c3o1 * (m_002 * m_120 + m_020 * m_102 + c4o1 * m_011 * m_111 + c2o1 * (m_101 * m_021 + m_110 * m_012)) +
            (m_120 + m_102)) * oneOverRho;
    m_212 = c_212 + c1o3 *
            (c3o1 * (m_002 * m_210 + m_200 * m_012 + c4o1 * m_101 * m_111 + c2o1 * (m_011 * m_201 + m_110 * m_102)) +
            (m_210 + m_012)) * oneOverRho;
    m_221 = c_221 + c1o3 *
            (c3o1 * (m_200 * m_021 + m_020 * m_201 + c4o1 * m_110 * m_111 + c2o1 * (m_101 * m_120 + m_011 * m_210)) +
            (m_021 + m_201)) * oneOverRho;

    //////////////////////////////////////////////////////////////////////////
    // 6.
    m_222 = c_222 - ((-c4o1 * m_111 * m_111 - (m_200 * m_022 + m_020 * m_202 + m_002 * m_220) -
                    c4o1 * (m_011 * m_211 + m_101 * m_121 + m_110 * m_112) -
                    c2o1 * (m_120 * m_102 + m_210 * m_012 + m_201 * m_021)) *
                    oneOverRho +
                    (c4o1 * (m_101 * m_101 * m_020 + m_011 * m_011 * m_200 + m_110 * m_110 * m_002) +
                    c2o1 * (m_200 * m_020 * m_002) + c16o1 * m_110 * m_101 * m_011) *
                    oneOverRho * oneOverRho -
                    c1o3 * (m_022 + m_202 + m_220) * oneOverRho - c1o9 * (m_200 + m_020 + m_002) * oneOverRho +
                    (c2o1 * (m_101 * m_101 + m_011 * m_011 + m_110 * m_110) +
                    (m_002 * m_020 + m_002 * m_200 + m_020 * m_200) + c1o3 * (m_002 + m_020 + m_200)) *
                    oneOverRho * oneOverRho * c2o3 +
                    c1o27 * ((drho * drho - drho) * oneOverRho * oneOverRho));

    ////////////////////////////////////////////////////////////////////////////////////
    //! -  Add acceleration (body force) to first order cumulants according to Eq. (85)-(87) in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    m_100 = -m_100;
    m_010 = -m_010;
    m_001 = -m_001;

    //Write to array here to distribute read/write
    if(writeMacroscopicVariables || turbulenceModel==TurbulenceModel::AMD)
    {
        rho[k_000] = drho;
        vx[k_000] = vvx;
        vy[k_000] = vvy;
        vz[k_000] = vvz;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardInverseChimeraWithK(m_000, m_100, m_200, vvx, vx2, c1o1, c1o1);
    backwardChimera(            m_010, m_110, m_210, vvx, vx2);
    backwardInverseChimeraWithK(m_020, m_120, m_220, vvx, vx2, c3o1, c1o3);
    backwardChimera(            m_001, m_101, m_201, vvx, vx2);
    backwardChimera(            m_011, m_111, m_211, vvx, vx2);
    backwardChimera(            m_021, m_121, m_221, vvx, vx2);
    backwardInverseChimeraWithK(m_002, m_102, m_202, vvx, vx2, c3o1, c1o3);
    backwardChimera(            m_012, m_112, m_212, vvx, vx2);
    backwardInverseChimeraWithK(m_022, m_122, m_222, vvx, vx2, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardInverseChimeraWithK(m_000, m_010, m_020, vvy, vy2, c6o1, c1o6);
    backwardChimera(            m_001, m_011, m_021, vvy, vy2);
    backwardInverseChimeraWithK(m_002, m_012, m_022, vvy, vy2, c18o1, c1o18);
    backwardInverseChimeraWithK(m_100, m_110, m_120, vvy, vy2, c3o2, c2o3);
    backwardChimera(            m_101, m_111, m_121, vvy, vy2);
    backwardInverseChimeraWithK(m_102, m_112, m_122, vvy, vy2, c9o2, c2o9);
    backwardInverseChimeraWithK(m_200, m_210, m_220, vvy, vy2, c6o1, c1o6);
    backwardChimera(            m_201, m_211, m_221, vvy, vy2);
    backwardInverseChimeraWithK(m_202, m_212, m_222, vvy, vy2, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardInverseChimeraWithK(m_000, m_001, m_002, vvz, vz2, c36o1, c1o36);
    backwardInverseChimeraWithK(m_010, m_011, m_012, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m_020, m_021, m_022, vvz, vz2, c36o1, c1o36);
    backwardInverseChimeraWithK(m_100, m_101, m_102, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m_110, m_111, m_112, vvz, vz2, c9o4, c4o9);
    backwardInverseChimeraWithK(m_120, m_121, m_122, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m_200, m_201, m_202, vvz, vz2, c36o1, c1o36);
    backwardInverseChimeraWithK(m_210, m_211, m_212, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m_220, m_221, m_222, vvz, vz2, c36o1, c1o36);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (dist.f[DIR_P00])[k_000] = f_M00;
    (dist.f[DIR_M00])[k_M00] = f_P00;
    (dist.f[DIR_0P0])[k_000] = f_0M0;
    (dist.f[DIR_0M0])[k_0M0] = f_0P0;
    (dist.f[DIR_00P])[k_000] = f_00M;
    (dist.f[DIR_00M])[k_00M] = f_00P;
    (dist.f[DIR_PP0])[k_000] = f_MM0;
    (dist.f[DIR_MM0])[k_MM0] = f_PP0;
    (dist.f[DIR_PM0])[k_0M0] = f_MP0;
    (dist.f[DIR_MP0])[k_M00] = f_PM0;
    (dist.f[DIR_P0P])[k_000] = f_M0M;
    (dist.f[DIR_M0M])[k_M0M] = f_P0P;
    (dist.f[DIR_P0M])[k_00M] = f_M0P;
    (dist.f[DIR_M0P])[k_M00] = f_P0M;
    (dist.f[DIR_0PP])[k_000] = f_0MM;
    (dist.f[DIR_0MM])[k_0MM] = f_0PP;
    (dist.f[DIR_0PM])[k_00M] = f_0MP;
    (dist.f[DIR_0MP])[k_0M0] = f_0PM;
    (dist.f[DIR_000])[k_000] = f_000;
    (dist.f[DIR_PPP])[k_000] = f_MMM;
    (dist.f[DIR_PMP])[k_0M0] = f_MPM;
    (dist.f[DIR_PPM])[k_00M] = f_MMP;
    (dist.f[DIR_PMM])[k_0MM] = f_MPP;
    (dist.f[DIR_MPP])[k_M00] = f_PMM;
    (dist.f[DIR_MMP])[k_MM0] = f_PPM;
    (dist.f[DIR_MPM])[k_M0M] = f_PMP;
    (dist.f[DIR_MMM])[k_MMM] = f_PPP;
}

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::AMD, true, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::Smagorinsky, true, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::QR, true, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::None, true, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::AMD, true, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::Smagorinsky, true, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::QR, true, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::None, true, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::AMD, false, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::Smagorinsky, false, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::QR, false, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::None, false, true > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::AMD, false, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::Smagorinsky, false, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::QR, false, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);

template __global__ void LB_Kernel_CumulantK17 < TurbulenceModel::None, false, false > ( real omega_in, uint* neighborX, uint* neighborY, uint* neighborZ, real* distributions, real* rho, real* vx, real* vy, real* vz, real* turbulentViscosity, real SGSconstant, unsigned long long numberOfLBnodes, int level, real* forces, real* bodyForceX, real* bodyForceY, real* bodyForceZ, real* quadricLimiters, bool isEvenTimestep, const uint *fluidNodeIndices, uint numberOfFluidNodes);
