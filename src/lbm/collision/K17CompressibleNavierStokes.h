
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
//! \brief Kernel for CumulantK17 including different turbulence models and options for local body forces and writing
//! macroscopic variables
//!
//! CumulantK17 kernel using chimera transformations and quartic limiters as present in Geier et al. (2017). Additional
//! options are three different eddy-viscosity turbulence models (Smagorinsky, AMD, QR) that can be set via the template
//! parameter turbulenceModel (with default TurbulenceModel::None). The kernel is executed separately for each subset of
//! fluid node indices with a different tag CollisionTemplate. For each subset, only the locally required options are
//! switched on ( \param writeMacroscopicVariables and/or \param applyBodyForce) in order to minimize memory accesses. The
//! default refers to the plain cumlant kernel (CollisionTemplate::Default). Nodes are added to subsets (taggedFluidNodes) in
//! Simulation::init using a corresponding tag with different values of CollisionTemplate. These subsets are provided by the
//! utilized PostCollisionInteractiors depending on they specific requirements (e.g. writeMacroscopicVariables for probes).

//=======================================================================================
#include <basics/constants/NumericConstants.h>

#include "lbm/constants/D3Q27.h"

#include "lbm/ChimeraTransformation.h"

#include "TurbulentViscosity.h"

#include "lbm/MacroscopicQuantities.h"

#include "CollisionParameter.h"

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#ifdef __CUDACC__
#define KERNEL_ABS abs
#else
#include <cmath>
#define KERNEL_ABS std::abs
#endif

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

//////////////////////////////////////////////////////////////////////////
//! Cumulant K17 Kernel is based on \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> and \ref <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017),
//! DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
////////////////////////////////////////////////////////////////////////////////
template <TurbulenceModel turbulenceModel>
__host__ __device__ void runK17CompressibleNavierStokes(CollisionParameter& parameter, MacroscopicValues& macroscopicValues, TurbulentViscosity& turbulentViscosity)
{
    auto& distribution = parameter.distribution;

    real& f000 = distribution[d000];
    real& fP00 = distribution[dP00];
    real& fM00 = distribution[dM00];
    real& f0P0 = distribution[d0P0];
    real& f0M0 = distribution[d0M0];
    real& f00P = distribution[d00P];
    real& f00M = distribution[d00M];
    real& fPP0 = distribution[dPP0];
    real& fMM0 = distribution[dMM0];
    real& fPM0 = distribution[dPM0];
    real& fMP0 = distribution[dMP0];
    real& fP0P = distribution[dP0P];
    real& fM0M = distribution[dM0M];
    real& fP0M = distribution[dP0M];
    real& fM0P = distribution[dM0P];
    real& f0PP = distribution[d0PP];
    real& f0MM = distribution[d0MM];
    real& f0PM = distribution[d0PM];
    real& f0MP = distribution[d0MP];
    real& fPPP = distribution[dPPP];
    real& fMPP = distribution[dMPP];
    real& fPMP = distribution[dPMP];
    real& fMMP = distribution[dMMP];
    real& fPPM = distribution[dPPM];
    real& fMPM = distribution[dMPM];
    real& fPMM = distribution[dPMM];
    real& fMMM = distribution[dMMM];

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Define aliases to use the same variable for the moments (m's):
    //!
    real& m111 = f000;
    real& m211 = fP00;
    real& m011 = fM00;
    real& m121 = f0P0;
    real& m101 = f0M0;
    real& m112 = f00P;
    real& m110 = f00M;
    real& m221 = fPP0;
    real& m001 = fMM0;
    real& m201 = fPM0;
    real& m021 = fMP0;
    real& m212 = fP0P;
    real& m010 = fM0M;
    real& m210 = fP0M;
    real& m012 = fM0P;
    real& m122 = f0PP;
    real& m100 = f0MM;
    real& m120 = f0PM;
    real& m102 = f0MP;
    real& m222 = fPPP;
    real& m022 = fMPP;
    real& m202 = fPMP;
    real& m002 = fMMP;
    real& m220 = fPPM;
    real& m020 = fMPM;
    real& m200 = fPMM;
    real& m000 = fMMM;

    //////////////////////////////////////////////////////(unsigned long)//////////////////////////////
    //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    real drho, oneOverRho, vvx, vvy, vvz;
    getCompressibleMacroscopicValues(distribution, drho, oneOverRho, vvx, vvy, vvz);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Add half of the acceleration (body force) to the velocity as in Eq. (42) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    vvx += parameter.forceX;
    vvy += parameter.forceY;
    vvz += parameter.forceZ;

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
    real quadricLimitP = parameter.quadricLimiter[0];
    real quadricLimitM = parameter.quadricLimiter[1];
    real quadricLimitD = parameter.quadricLimiter[2];
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from well conditioned distributions to central moments as defined in Appendix J in \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (6)-(14) in \ref <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    forwardInverseChimeraWithK(fMMM, fMM0, fMMP, vvz, vz2, c36o1, c1o36);
    forwardInverseChimeraWithK(fM0M, fM00, fM0P, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(fMPM, fMP0, fMPP, vvz, vz2, c36o1, c1o36);
    forwardInverseChimeraWithK(f0MM, f0M0, f0MP, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(f00M, f000, f00P, vvz, vz2, c9o4,  c4o9);
    forwardInverseChimeraWithK(f0PM, f0P0, f0PP, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(fPMM, fPM0, fPMP, vvz, vz2, c36o1, c1o36);
    forwardInverseChimeraWithK(fP0M, fP00, fP0P, vvz, vz2, c9o1,  c1o9);
    forwardInverseChimeraWithK(fPPM, fPP0, fPPP, vvz, vz2, c36o1, c1o36);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    forwardInverseChimeraWithK(fMMM, fM0M, fMPM, vvy, vy2, c6o1,  c1o6);
    forwardChimera(            fMM0, fM00, fMP0, vvy, vy2);
    forwardInverseChimeraWithK(fMMP, fM0P, fMPP, vvy, vy2, c18o1, c1o18);
    forwardInverseChimeraWithK(f0MM, f00M, f0PM, vvy, vy2, c3o2,  c2o3);
    forwardChimera(            f0M0, f000, f0P0, vvy, vy2);
    forwardInverseChimeraWithK(f0MP, f00P, f0PP, vvy, vy2, c9o2,  c2o9);
    forwardInverseChimeraWithK(fPMM, fP0M, fPPM, vvy, vy2, c6o1,  c1o6);
    forwardChimera(            fPM0, fP00, fPP0, vvy, vy2);
    forwardInverseChimeraWithK(fPMP, fP0P, fPPP, vvy, vy2, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    forwardInverseChimeraWithK(fMMM, f0MM, fPMM, vvx, vx2, c1o1, c1o1);
    forwardChimera(            fM0M, f00M, fP0M, vvx, vx2);
    forwardInverseChimeraWithK(fMPM, f0PM, fPPM, vvx, vx2, c3o1, c1o3);
    forwardChimera(            fMM0, f0M0, fPM0, vvx, vx2);
    forwardChimera(            fM00, f000, fP00, vvx, vx2);
    forwardChimera(            fMP0, f0P0, fPP0, vvx, vx2);
    forwardInverseChimeraWithK(fMMP, f0MP, fPMP, vvx, vx2, c3o1, c1o3);
    forwardChimera(            fM0P, f00P, fP0P, vvx, vx2);
    forwardInverseChimeraWithK(fMPP, f0PP, fPPP, vvx, vx2, c3o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Setting relaxation rates for non-hydrodynamic cumulants (default values). Variable names and equations
    //! according to <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!  => [NAME IN PAPER]=[NAME IN CODE]=[DEFAULT VALUE].
    //!  - Trace of second order cumulants \f$ C{200}+C{020}+C{002} \f$ used to adjust bulk
    //!  viscosity:\f$\omega2=OxxPyyPzz=1.0 \f$.
    //!  - Third order cumulants \f$ C{120}+C{102}, C{210}+C{012}, C{201}+C{021} \f$: \f$ \omega3=OxyyPxzz
    //!  \f$ set according to Eq. (111) with simplifications assuming \f$ \omega2=1.0\f$.
    //!  - Third order cumulants \f$ C{120}-C{102}, C{210}-C{012}, C{201}-C{021} \f$: \f$ \omega4 = OxyyMxzz
    //!  \f$ set according to Eq. (112) with simplifications assuming \f$ \omega2 = 1.0\f$.
    //!  - Third order cumulants \f$ C{111} \f$: \f$ \omega5 = Oxyz \f$ set according to Eq. (113) with
    //!  simplifications assuming \f$ \omega2 = 1.0\f$  (modify for different bulk viscosity).
    //!  - Fourth order cumulants \f$ C{220}, C{202}, C{022}, C{211}, C{121}, C{112} \f$: for simplification
    //!  all set to the same default value \f$ \omega6=\omega7=\omega8=O4=1.0 \f$.
    //!  - Fifth order cumulants \f$ C{221}, C{212}, C{122}\f$: \f$\omega9=O5=1.0\f$.
    //!  - Sixth order cumulant \f$ C{222}\f$: \f$\omega{10}=O6=1.0\f$.
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Calculate modified omega with turbulent viscosity
    //!
    const real omega = turbulenceModel == TurbulenceModel::None ? parameter.omega : vf::lbm::calculateOmegaWithturbulentViscosity(parameter.omega, turbulentViscosity.value);
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
    //! - A and d00M: parameters for fourth order convergence of the diffusion term according to Eq. (115) and (116)
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
    real cm211 = m211 - ((m200 + c1o3) * m011 + c2o1 * m110 * m101) * oneOverRho;
    real cm121 = m121 - ((m020 + c1o3) * m101 + c2o1 * m110 * m011) * oneOverRho;
    real cm112 = m112 - ((m002 + c1o3) * m110 + c2o1 * m101 * m011) * oneOverRho;

    real cm220 = m220 - (((m200 * m020 + c2o1 * m110 * m110) + c1o3 * (m200 + m020)) * oneOverRho - c1o9 * (drho * oneOverRho));
    real cm202 = m202 - (((m200 * m002 + c2o1 * m101 * m101) + c1o3 * (m200 + m002)) * oneOverRho - c1o9 * (drho * oneOverRho));
    real cm022 = m022 - (((m002 * m020 + c2o1 * m011 * m011) + c1o3 * (m002 + m020)) * oneOverRho - c1o9 * (drho * oneOverRho));
    ////////////////////////////////////////////////////////////
    // 5.
    real cm122 =
        m122 - ((m002 * m120 + m020 * m102 + c4o1 * m011 * m111 + c2o1 * (m101 * m021 + m110 * m012)) +
                c1o3 * (m120 + m102)) *
                oneOverRho;
    real cm212 =
        m212 - ((m002 * m210 + m200 * m012 + c4o1 * m101 * m111 + c2o1 * (m011 * m201 + m110 * m102)) +
                c1o3 * (m210 + m012)) *
                oneOverRho;
    real cm221 =
        m221 - ((m200 * m021 + m020 * m201 + c4o1 * m110 * m111 + c2o1 * (m101 * m120 + m011 * m210)) +
                c1o3 * (m021 + m201)) *
                oneOverRho;
    ////////////////////////////////////////////////////////////
    // 6.
    real cm222 = m222 + ((-c4o1 * m111 * m111 - (m200 * m022 + m020 * m202 + m002 * m220) -
                            c4o1 * (m011 * m211 + m101 * m121 + m110 * m112) -
                            c2o1 * (m120 * m102 + m210 * m012 + m201 * m021)) *
                            oneOverRho +
                        (c4o1 * (m101 * m101 * m020 + m011 * m011 * m200 + m110 * m110 * m002) +
                            c2o1 * (m200 * m020 * m002) + c16o1 * m110 * m101 * m011) *
                            oneOverRho * oneOverRho -
                            c1o3 * (m022 + m202 + m220) * oneOverRho - c1o9 * (m200 + m020 + m002) * oneOverRho +
                        (c2o1 * (m101 * m101 + m011 * m011 + m110 * m110) +
                            (m002 * m020 + m002 * m200 + m020 * m200) + c1o3 * (m002 + m020 + m200)) *
                            oneOverRho * oneOverRho * c2o3 +
                            c1o27 * ((drho * drho - drho) * oneOverRho * oneOverRho));

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute linear combinations of second and third order cumulants
    //!
    ////////////////////////////////////////////////////////////
    // 2.
    real mxxPyyPzz = m200 + m020 + m002;
    real mxxMyy    = m200 - m020;
    real mxxMzz    = m200 - m002;
    ////////////////////////////////////////////////////////////
    // 3.
    real mxxyPyzz = m210 + m012;
    real mxxyMyzz = m210 - m012;

    real mxxzPyyz = m201 + m021;
    real mxxzMyyz = m201 - m021;

    real mxyyPxzz = m120 + m102;
    real mxyyMxzz = m120 - m102;

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
    real Dxy = -c3o1 * omega * m110;
    real Dxz = -c3o1 * omega * m101;
    real Dyz = -c3o1 * omega * m011;
    real dxux = c1o2 * (-omega) * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (m000 - mxxPyyPzz);
    real dyuy = dxux + omega * c3o2 * mxxMyy;
    real dzuz = dxux + omega * c3o2 * mxxMzz;

    ////////////////////////////////////////////////////////////////////////////////////
    switch (turbulenceModel) {
        case TurbulenceModel::None:
        case TurbulenceModel::AMD: // AMD is computed in separate kernel
            break;
        case TurbulenceModel::Smagorinsky:
            turbulentViscosity.value =
                calcTurbulentViscositySmagorinsky(turbulentViscosity.SGSconstant, dxux, dyuy, dzuz, Dxy, Dxz, Dyz);
            break;
        case TurbulenceModel::QR:
            turbulentViscosity.value =
                calcTurbulentViscosityQR(turbulentViscosity.SGSconstant, dxux, dyuy, dzuz, Dxy, Dxz, Dyz);
            break;
        default:
            break;
    }
    ////////////////////////////////////////////////////////////
    //! - Relaxation of second order cumulants with correction terms according to Eq. (33)-(35) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!
    mxxPyyPzz += OxxPyyPzz * (m000 - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
    mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
    mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

    ////////////////////////////////////////////////////////////////////////////////////
    ////no correction
    // mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);
    // mxxMyy += -(-omega) * (-mxxMyy);
    // mxxMzz += -(-omega) * (-mxxMzz);
    //////////////////////////////////////////////////////////////////////////
    m011 += omega * (-m011);
    m101 += omega * (-m101);
    m110 += omega * (-m110);

    ////////////////////////////////////////////////////////////////////////////////////
    // relax
    //////////////////////////////////////////////////////////////////////////
    // incl. limiter
    //! - Relaxation of third order cumulants including limiter according to Eq. (116)-(123)
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!
    real wadjust = Oxyz + (c1o1 - Oxyz) * KERNEL_ABS(m111) / (KERNEL_ABS(m111) + quadricLimitD);
    m111 += wadjust * (-m111);
    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * KERNEL_ABS(mxxyPyzz) / (KERNEL_ABS(mxxyPyzz) + quadricLimitP);
    mxxyPyzz += wadjust * (-mxxyPyzz);
    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * KERNEL_ABS(mxxyMyzz) / (KERNEL_ABS(mxxyMyzz) + quadricLimitM);
    mxxyMyzz += wadjust * (-mxxyMyzz);
    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * KERNEL_ABS(mxxzPyyz) / (KERNEL_ABS(mxxzPyyz) + quadricLimitP);
    mxxzPyyz += wadjust * (-mxxzPyyz);
    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * KERNEL_ABS(mxxzMyyz) / (KERNEL_ABS(mxxzMyyz) + quadricLimitM);
    mxxzMyyz += wadjust * (-mxxzMyyz);
    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * KERNEL_ABS(mxyyPxzz) / (KERNEL_ABS(mxyyPxzz) + quadricLimitP);
    mxyyPxzz += wadjust * (-mxyyPxzz);
    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * KERNEL_ABS(mxyyMxzz) / (KERNEL_ABS(mxyyMxzz) + quadricLimitM);
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
    m200 = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
    m020 = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz);
    m002 = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

    m210 = ( mxxyMyzz + mxxyPyzz) * c1o2;
    m012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
    m201 = ( mxxzMyyz + mxxzPyyz) * c1o2;
    m021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
    m120 = ( mxyyMxzz + mxyyPxzz) * c1o2;
    m102 = (-mxyyMxzz + mxyyPxzz) * c1o2;
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // 4.
    // no limiter
    //! - Relax fourth order cumulants to modified equilibrium for fourth order convergence of diffusion according
    //! to Eq. (43)-(48) <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!
    cm022 = -O4 * (c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * factorA + (c1o1 - O4) * (cm022);
    cm202 = -O4 * (c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * factorA + (c1o1 - O4) * (cm202);
    cm220 = -O4 * (c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * factorA + (c1o1 - O4) * (cm220);
    cm112 = -O4 * (c1o1 / omega - c1o2) * Dxy           * c1o3 * factorB + (c1o1 - O4) * (cm112);
    cm121 = -O4 * (c1o1 / omega - c1o2) * Dxz           * c1o3 * factorB + (c1o1 - O4) * (cm121);
    cm211 = -O4 * (c1o1 / omega - c1o2) * Dyz           * c1o3 * factorB + (c1o1 - O4) * (cm211);


    //////////////////////////////////////////////////////////////////////////
    // 5.
    cm122 += O5 * (-cm122);
    cm212 += O5 * (-cm212);
    cm221 += O5 * (-cm221);

    //////////////////////////////////////////////////////////////////////////
    // 6.
    cm222 += O6 * (-cm222);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute central moments from post collision cumulants according to Eq. (53)-(56) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //!

    //////////////////////////////////////////////////////////////////////////
    // 4.
    m211 = cm211 + c1o3 * ((c3o1 * m200 + c1o1) * m011 + c6o1 * m110 * m101) * oneOverRho;
    m121 = cm121 + c1o3 * ((c3o1 * m020 + c1o1) * m101 + c6o1 * m110 * m011) * oneOverRho;
    m112 = cm112 + c1o3 * ((c3o1 * m002 + c1o1) * m110 + c6o1 * m101 * m011) * oneOverRho;

    m220 =
        cm220 + (((m200 * m020 + c2o1 * m110 * m110) * c9o1 + c3o1 * (m200 + m020)) * oneOverRho - (drho * oneOverRho)) * c1o9;
    m202 =
        cm202 + (((m200 * m002 + c2o1 * m101 * m101) * c9o1 + c3o1 * (m200 + m002)) * oneOverRho - (drho * oneOverRho)) * c1o9;
    m022 =
        cm022 + (((m002 * m020 + c2o1 * m011 * m011) * c9o1 + c3o1 * (m002 + m020)) * oneOverRho - (drho * oneOverRho)) * c1o9;

    //////////////////////////////////////////////////////////////////////////
    // 5.
    m122 = cm122 + c1o3 *
            (c3o1 * (m002 * m120 + m020 * m102 + c4o1 * m011 * m111 + c2o1 * (m101 * m021 + m110 * m012)) +
            (m120 + m102)) * oneOverRho;
    m212 = cm212 + c1o3 *
            (c3o1 * (m002 * m210 + m200 * m012 + c4o1 * m101 * m111 + c2o1 * (m011 * m201 + m110 * m102)) +
            (m210 + m012)) * oneOverRho;
    m221 = cm221 + c1o3 *
            (c3o1 * (m200 * m021 + m020 * m201 + c4o1 * m110 * m111 + c2o1 * (m101 * m120 + m011 * m210)) +
            (m021 + m201)) * oneOverRho;

    //////////////////////////////////////////////////////////////////////////
    // 6.
    m222 = cm222 - ((-c4o1 * m111 * m111 - (m200 * m022 + m020 * m202 + m002 * m220) -
                    c4o1 * (m011 * m211 + m101 * m121 + m110 * m112) -
                    c2o1 * (m120 * m102 + m210 * m012 + m201 * m021)) *
                    oneOverRho +
                    (c4o1 * (m101 * m101 * m020 + m011 * m011 * m200 + m110 * m110 * m002) +
                    c2o1 * (m200 * m020 * m002) + c16o1 * m110 * m101 * m011) *
                    oneOverRho * oneOverRho -
                    c1o3 * (m022 + m202 + m220) * oneOverRho - c1o9 * (m200 + m020 + m002) * oneOverRho +
                    (c2o1 * (m101 * m101 + m011 * m011 + m110 * m110) +
                    (m002 * m020 + m002 * m200 + m020 * m200) + c1o3 * (m002 + m020 + m200)) *
                    oneOverRho * oneOverRho * c2o3 +
                    c1o27 * ((drho * drho - drho) * oneOverRho * oneOverRho));

    ////////////////////////////////////////////////////////////////////////////////////
    //! -  Add acceleration (body force) to first order cumulants according to Eq. (85)-(87) in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
    //!
    m100 = -m100;
    m010 = -m010;
    m001 = -m001;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardInverseChimeraWithK(m000, m100, m200, vvx, vx2, c1o1, c1o1);
    backwardChimera(            m010, m110, m210, vvx, vx2);
    backwardInverseChimeraWithK(m020, m120, m220, vvx, vx2, c3o1, c1o3);
    backwardChimera(            m001, m101, m201, vvx, vx2);
    backwardChimera(            m011, m111, m211, vvx, vx2);
    backwardChimera(            m021, m121, m221, vvx, vx2);
    backwardInverseChimeraWithK(m002, m102, m202, vvx, vx2, c3o1, c1o3);
    backwardChimera(            m012, m112, m212, vvx, vx2);
    backwardInverseChimeraWithK(m022, m122, m222, vvx, vx2, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardInverseChimeraWithK(m000, m010, m020, vvy, vy2, c6o1, c1o6);
    backwardChimera(            m001, m011, m021, vvy, vy2);
    backwardInverseChimeraWithK(m002, m012, m022, vvy, vy2, c18o1, c1o18);
    backwardInverseChimeraWithK(m100, m110, m120, vvy, vy2, c3o2, c2o3);
    backwardChimera(            m101, m111, m121, vvy, vy2);
    backwardInverseChimeraWithK(m102, m112, m122, vvy, vy2, c9o2, c2o9);
    backwardInverseChimeraWithK(m200, m210, m220, vvy, vy2, c6o1, c1o6);
    backwardChimera(            m201, m211, m221, vvy, vy2);
    backwardInverseChimeraWithK(m202, m212, m222, vvy, vy2, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardInverseChimeraWithK(m000, m001, m002, vvz, vz2, c36o1, c1o36);
    backwardInverseChimeraWithK(m010, m011, m012, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m020, m021, m022, vvz, vz2, c36o1, c1o36);
    backwardInverseChimeraWithK(m100, m101, m102, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m110, m111, m112, vvz, vz2, c9o4, c4o9);
    backwardInverseChimeraWithK(m120, m121, m122, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m200, m201, m202, vvz, vz2, c36o1, c1o36);
    backwardInverseChimeraWithK(m210, m211, m212, vvz, vz2, c9o1, c1o9);
    backwardInverseChimeraWithK(m220, m221, m222, vvz, vz2, c36o1, c1o36);

    distribution[dP00] = fP00;
    distribution[dM00] = fM00;
    distribution[d0P0] = f0P0;
    distribution[d0M0] = f0M0;
    distribution[d00P] = f00P;
    distribution[d00M] = f00M;
    distribution[dPP0] = fPP0;
    distribution[dMM0] = fMM0;
    distribution[dPM0] = fPM0;
    distribution[dMP0] = fMP0;
    distribution[dP0P] = fP0P;
    distribution[dM0M] = fM0M;
    distribution[dP0M] = fP0M;
    distribution[dM0P] = fM0P;
    distribution[d0PP] = f0PP;
    distribution[d0MM] = f0MM;
    distribution[d0PM] = f0PM;
    distribution[d0MP] = f0MP;
    distribution[d000] = f000;
    distribution[dPPP] = fPPP;
    distribution[dPMP] = fPMP;
    distribution[dPPM] = fPPM;
    distribution[dPMM] = fPMM;
    distribution[dMPP] = fMPP;
    distribution[dMMP] = fMMP;
    distribution[dMPM] = fMPM;
    distribution[dMMM] = fMMM;

    macroscopicValues.rho = drho;
    macroscopicValues.vx = vvx;
    macroscopicValues.vy = vvy;
    macroscopicValues.vz = vvz;
}

} // namespace vf::lbm
