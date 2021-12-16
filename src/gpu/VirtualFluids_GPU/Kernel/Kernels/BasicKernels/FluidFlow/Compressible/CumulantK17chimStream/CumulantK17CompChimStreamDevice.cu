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
//! \file Cumulant27chim.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
#include "Kernel/ChimeraTransformation.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_Kernel_CumulantK17CompChimStream(
	real omega,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	unsigned long size_Mat,
	int level,
	real* forces,
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
    const unsigned x = threadIdx.x;
    const unsigned y = blockIdx.x;
    const unsigned z = blockIdx.y;

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k_thread = nx * (ny * z + y) + x;

    //////////////////////////////////////////////////////////////////////////
    // run for all indices in fluidNodeIndices
    if (k_thread < numberOfFluidNodes) {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>

        const unsigned k = fluidNodeIndices[k_thread];

        Distributions27 dist;
        if (isEvenTimestep) {
            dist.f[dirE]    = &distributions[dirE * size_Mat];
            dist.f[dirW]    = &distributions[dirW * size_Mat];
            dist.f[dirN]    = &distributions[dirN * size_Mat];
            dist.f[dirS]    = &distributions[dirS * size_Mat];
            dist.f[dirT]    = &distributions[dirT * size_Mat];
            dist.f[dirB]    = &distributions[dirB * size_Mat];
            dist.f[dirNE]   = &distributions[dirNE * size_Mat];
            dist.f[dirSW]   = &distributions[dirSW * size_Mat];
            dist.f[dirSE]   = &distributions[dirSE * size_Mat];
            dist.f[dirNW]   = &distributions[dirNW * size_Mat];
            dist.f[dirTE]   = &distributions[dirTE * size_Mat];
            dist.f[dirBW]   = &distributions[dirBW * size_Mat];
            dist.f[dirBE]   = &distributions[dirBE * size_Mat];
            dist.f[dirTW]   = &distributions[dirTW * size_Mat];
            dist.f[dirTN]   = &distributions[dirTN * size_Mat];
            dist.f[dirBS]   = &distributions[dirBS * size_Mat];
            dist.f[dirBN]   = &distributions[dirBN * size_Mat];
            dist.f[dirTS]   = &distributions[dirTS * size_Mat];
            dist.f[dirZERO] = &distributions[dirZERO * size_Mat];
            dist.f[dirTNE]  = &distributions[dirTNE * size_Mat];
            dist.f[dirTSW]  = &distributions[dirTSW * size_Mat];
            dist.f[dirTSE]  = &distributions[dirTSE * size_Mat];
            dist.f[dirTNW]  = &distributions[dirTNW * size_Mat];
            dist.f[dirBNE]  = &distributions[dirBNE * size_Mat];
            dist.f[dirBSW]  = &distributions[dirBSW * size_Mat];
            dist.f[dirBSE]  = &distributions[dirBSE * size_Mat];
            dist.f[dirBNW]  = &distributions[dirBNW * size_Mat];
        } else {
            dist.f[dirW]    = &distributions[dirE * size_Mat];
            dist.f[dirE]    = &distributions[dirW * size_Mat];
            dist.f[dirS]    = &distributions[dirN * size_Mat];
            dist.f[dirN]    = &distributions[dirS * size_Mat];
            dist.f[dirB]    = &distributions[dirT * size_Mat];
            dist.f[dirT]    = &distributions[dirB * size_Mat];
            dist.f[dirSW]   = &distributions[dirNE * size_Mat];
            dist.f[dirNE]   = &distributions[dirSW * size_Mat];
            dist.f[dirNW]   = &distributions[dirSE * size_Mat];
            dist.f[dirSE]   = &distributions[dirNW * size_Mat];
            dist.f[dirBW]   = &distributions[dirTE * size_Mat];
            dist.f[dirTE]   = &distributions[dirBW * size_Mat];
            dist.f[dirTW]   = &distributions[dirBE * size_Mat];
            dist.f[dirBE]   = &distributions[dirTW * size_Mat];
            dist.f[dirBS]   = &distributions[dirTN * size_Mat];
            dist.f[dirTN]   = &distributions[dirBS * size_Mat];
            dist.f[dirTS]   = &distributions[dirBN * size_Mat];
            dist.f[dirBN]   = &distributions[dirTS * size_Mat];
            dist.f[dirZERO] = &distributions[dirZERO * size_Mat];
            dist.f[dirBSW]  = &distributions[dirTNE * size_Mat];
            dist.f[dirBNE]  = &distributions[dirTSW * size_Mat];
            dist.f[dirBNW]  = &distributions[dirTSE * size_Mat];
            dist.f[dirBSE]  = &distributions[dirTNW * size_Mat];
            dist.f[dirTSW]  = &distributions[dirBNE * size_Mat];
            dist.f[dirTNE]  = &distributions[dirBSW * size_Mat];
            dist.f[dirTNW]  = &distributions[dirBSE * size_Mat];
            dist.f[dirTSE]  = &distributions[dirBNW * size_Mat];
        }
        ////////////////////////////////////////////////////////////////////////////////
        //! - Set neighbor indices (necessary for indirect addressing)
        uint kw   = neighborX[k];
        uint ks   = neighborY[k];
        uint kb   = neighborZ[k];
        uint ksw  = neighborY[kw];
        uint kbw  = neighborZ[kw];
        uint kbs  = neighborZ[ks];
        uint kbsw = neighborZ[ksw];
        ////////////////////////////////////////////////////////////////////////////////////
        //! - Set local distributions
        //!
        real mfcbb = (dist.f[dirE])[k];
        real mfabb = (dist.f[dirW])[kw];
        real mfbcb = (dist.f[dirN])[k];
        real mfbab = (dist.f[dirS])[ks];
        real mfbbc = (dist.f[dirT])[k];
        real mfbba = (dist.f[dirB])[kb];
        real mfccb = (dist.f[dirNE])[k];
        real mfaab = (dist.f[dirSW])[ksw];
        real mfcab = (dist.f[dirSE])[ks];
        real mfacb = (dist.f[dirNW])[kw];
        real mfcbc = (dist.f[dirTE])[k];
        real mfaba = (dist.f[dirBW])[kbw];
        real mfcba = (dist.f[dirBE])[kb];
        real mfabc = (dist.f[dirTW])[kw];
        real mfbcc = (dist.f[dirTN])[k];
        real mfbaa = (dist.f[dirBS])[kbs];
        real mfbca = (dist.f[dirBN])[kb];
        real mfbac = (dist.f[dirTS])[ks];
        real mfbbb = (dist.f[dirZERO])[k];
        real mfccc = (dist.f[dirTNE])[k];
        real mfaac = (dist.f[dirTSW])[ksw];
        real mfcac = (dist.f[dirTSE])[ks];
        real mfacc = (dist.f[dirTNW])[kw];
        real mfcca = (dist.f[dirBNE])[kb];
        real mfaaa = (dist.f[dirBSW])[kbsw];
        real mfcaa = (dist.f[dirBSE])[kbs];
        real mfaca = (dist.f[dirBNW])[kbw];
        ////////////////////////////////////////////////////////////////////////////////////
        //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
        //!
        real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                     (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) +
                      ((mfacb + mfcab) + (mfaab + mfccb))) +
                     ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) +
                    mfbbb;

        real rho   = c1o1 + drho;
        real OOrho = c1o1 / rho;

        real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                    (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) + (mfcbb - mfabb)) *
                   OOrho;
        real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                    (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) + (mfbcb - mfbab)) *
                   OOrho;
        real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                    (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) + (mfbbc - mfbba)) *
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
        forwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36);
        forwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36);
        forwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c9o4, c4o9);
        forwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36);
        forwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c9o1, c1o9);
        forwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36);

        ////////////////////////////////////////////////////////////////////////////////////
        // Y - Dir
        forwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c6o1, c1o6);
        forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
        forwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18);
        forwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c3o2, c2o3);
        forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
        forwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c9o2, c2o9);
        forwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c6o1, c1o6);
        forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
        forwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18);

        ////////////////////////////////////////////////////////////////////////////////////
        // X - Dir
        forwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
        forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
        forwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3);
        forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
        forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
        forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
        forwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3);
        forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
        forwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c3o1, c1o9);

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
        //! - A and B: parameters for fourth order convergence of the diffusion term according to Eq. (114) and (115)
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a> with simplifications assuming \f$ \omega_2 = 1.0 \f$ (modify for
        //! different bulk viscosity).
        //!
        real A = (c4o1 + c2o1 * omega - c3o1 * omega * omega) / (c2o1 - c7o1 * omega + c5o1 * omega * omega);
        real B = (c4o1 + c28o1 * omega - c14o1 * omega * omega) / (c6o1 - c21o1 * omega + c15o1 * omega * omega);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Compute cumulants from central moments according to Eq. (20)-(23) in
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        ////////////////////////////////////////////////////////////
        // 4.
        real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) * OOrho;
        real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) * OOrho;
        real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) * OOrho;

        real CUMcca =
            mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) * OOrho - c1o9 * (drho * OOrho));
        real CUMcac =
            mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) * OOrho - c1o9 * (drho * OOrho));
        real CUMacc =
            mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) * OOrho - c1o9 * (drho * OOrho));
        ////////////////////////////////////////////////////////////
        // 5.
        real CUMbcc =
            mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) +
                     c1o3 * (mfbca + mfbac)) *
                        OOrho;
        real CUMcbc =
            mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) +
                     c1o3 * (mfcba + mfabc)) *
                        OOrho;
        real CUMccb =
            mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) +
                     c1o3 * (mfacb + mfcab)) *
                        OOrho;
        ////////////////////////////////////////////////////////////
        // 6.
        real CUMccc = mfccc + ((-c4o1 * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
                                c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc) -
                                c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) *
                                   OOrho +
                               (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
                                c2o1 * (mfcaa * mfaca * mfaac) + c16o1 * mfbba * mfbab * mfabb) *
                                   OOrho * OOrho -
                               c1o3 * (mfacc + mfcac + mfcca) * OOrho - c1o9 * (mfcaa + mfaca + mfaac) * OOrho +
                               (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
                                (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 * (mfaac + mfaca + mfcaa)) *
                                   OOrho * OOrho * c2o3 +
                               c1o27 * ((drho * drho - drho) * OOrho * OOrho));

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Compute linear combinations of second and third order cumulants
        //!
        ////////////////////////////////////////////////////////////
        // 2.
        real mxxPyyPzz = mfcaa + mfaca + mfaac;
        real mxxMyy    = mfcaa - mfaca;
        real mxxMzz    = mfcaa - mfaac;
        ////////////////////////////////////////////////////////////
        // 3.
        real mxxyPyzz = mfcba + mfabc;
        real mxxyMyzz = mfcba - mfabc;

        real mxxzPyyz = mfcab + mfacb;
        real mxxzMyyz = mfcab - mfacb;

        real mxyyPxzz = mfbca + mfbac;
        real mxyyMxzz = mfbca - mfbac;

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
        real Dxy  = -c3o1 * omega * mfbba;
        real Dxz  = -c3o1 * omega * mfbab;
        real Dyz  = -c3o1 * omega * mfabb;
        real dxux = c1o2 * (-omega) * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
        real dyuy = dxux + omega * c3o2 * mxxMyy;
        real dzuz = dxux + omega * c3o2 * mxxMzz;
        ////////////////////////////////////////////////////////////
        //! - Relaxation of second order cumulants with correction terms according to Eq. (33)-(35) in
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        mxxPyyPzz +=
            OxxPyyPzz * (mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
        mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
        mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

        ////////////////////////////////////////////////////////////////////////////////////
        ////no correction
        // mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);
        // mxxMyy += -(-omega) * (-mxxMyy);
        // mxxMzz += -(-omega) * (-mxxMzz);
        //////////////////////////////////////////////////////////////////////////
        mfabb += omega * (-mfabb);
        mfbab += omega * (-mfbab);
        mfbba += omega * (-mfbba);

        ////////////////////////////////////////////////////////////////////////////////////
        // relax
        //////////////////////////////////////////////////////////////////////////
        // incl. limiter
        //! - Relaxation of third order cumulants including limiter according to Eq. (116)-(123)
        //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        wadjust = Oxyz + (c1o1 - Oxyz) * abs(mfbbb) / (abs(mfbbb) + qudricLimitD);
        mfbbb += wadjust * (-mfbbb);
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
        mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
        mfaca = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz);
        mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

        mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
        mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
        mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
        mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
        mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
        mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
        //////////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////////
        // 4.
        // no limiter
        //! - Relax fourth order cumulants to modified equilibrium for fourth order convergence of diffusion according
        //! to Eq. (43)-(48) <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
        //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
        //!
        CUMacc = -O4 * (c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMacc);
        CUMcac = -O4 * (c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMcac);
        CUMcca = -O4 * (c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (c1o1 - O4) * (CUMcca);
        CUMbbc = -O4 * (c1o1 / omega - c1o2) * Dxy * c1o3 * B + (c1o1 - O4) * (CUMbbc);
        CUMbcb = -O4 * (c1o1 / omega - c1o2) * Dxz * c1o3 * B + (c1o1 - O4) * (CUMbcb);
        CUMcbb = -O4 * (c1o1 / omega - c1o2) * Dyz * c1o3 * B + (c1o1 - O4) * (CUMcbb);

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
        mfcbb = CUMcbb + c1o3 * ((c3o1 * mfcaa + c1o1) * mfabb + c6o1 * mfbba * mfbab) * OOrho;
        mfbcb = CUMbcb + c1o3 * ((c3o1 * mfaca + c1o1) * mfbab + c6o1 * mfbba * mfabb) * OOrho;
        mfbbc = CUMbbc + c1o3 * ((c3o1 * mfaac + c1o1) * mfbba + c6o1 * mfbab * mfabb) * OOrho;

        mfcca =
            CUMcca +
            (((mfcaa * mfaca + c2o1 * mfbba * mfbba) * c9o1 + c3o1 * (mfcaa + mfaca)) * OOrho - (drho * OOrho)) * c1o9;
        mfcac =
            CUMcac +
            (((mfcaa * mfaac + c2o1 * mfbab * mfbab) * c9o1 + c3o1 * (mfcaa + mfaac)) * OOrho - (drho * OOrho)) * c1o9;
        mfacc =
            CUMacc +
            (((mfaac * mfaca + c2o1 * mfabb * mfabb) * c9o1 + c3o1 * (mfaac + mfaca)) * OOrho - (drho * OOrho)) * c1o9;

        //////////////////////////////////////////////////////////////////////////
        // 5.
        mfbcc = CUMbcc + c1o3 *
                             (c3o1 * (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb +
                                      c2o1 * (mfbab * mfacb + mfbba * mfabc)) +
                              (mfbca + mfbac)) *
                             OOrho;
        mfcbc = CUMcbc + c1o3 *
                             (c3o1 * (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb +
                                      c2o1 * (mfabb * mfcab + mfbba * mfbac)) +
                              (mfcba + mfabc)) *
                             OOrho;
        mfccb = CUMccb + c1o3 *
                             (c3o1 * (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb +
                                      c2o1 * (mfbab * mfbca + mfabb * mfcba)) +
                              (mfacb + mfcab)) *
                             OOrho;

        //////////////////////////////////////////////////////////////////////////
        // 6.
        mfccc = CUMccc - ((-c4o1 * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
                           c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc) -
                           c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) *
                              OOrho +
                          (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
                           c2o1 * (mfcaa * mfaca * mfaac) + c16o1 * mfbba * mfbab * mfabb) *
                              OOrho * OOrho -
                          c1o3 * (mfacc + mfcac + mfcca) * OOrho - c1o9 * (mfcaa + mfaca + mfaac) * OOrho +
                          (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
                           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 * (mfaac + mfaca + mfcaa)) *
                              OOrho * OOrho * c2o3 +
                          c1o27 * ((drho * drho - drho) * OOrho * OOrho));

        ////////////////////////////////////////////////////////////////////////////////////
        //! -  Add acceleration (body force) to first order cumulants according to Eq. (85)-(87) in
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
        //!
        mfbaa = -mfbaa;
        mfaba = -mfaba;
        mfaab = -mfaab;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
        //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
        //! ]</b></a>
        //!
        ////////////////////////////////////////////////////////////////////////////////////
        // X - Dir
        backwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
        backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
        backwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3);
        backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
        backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
        backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
        backwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3);
        backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
        backwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c9o1, c1o9);

        ////////////////////////////////////////////////////////////////////////////////////
        // Y - Dir
        backwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c6o1, c1o6);
        backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
        backwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18);
        backwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c3o2, c2o3);
        backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
        backwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c9o2, c2o9);
        backwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c6o1, c1o6);
        backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
        backwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18);

        ////////////////////////////////////////////////////////////////////////////////////
        // Z - Dir
        backwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c9o4, c4o9);
        backwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c9o1, c1o9);
        backwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Write distributions: style of reading and writing the distributions from/to
        //! stored arrays dependent on timestep is based on the esoteric twist algorithm
        //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        (dist.f[dirE])[k]      = mfabb;
        (dist.f[dirW])[kw]     = mfcbb;
        (dist.f[dirN])[k]      = mfbab;
        (dist.f[dirS])[ks]     = mfbcb;
        (dist.f[dirT])[k]      = mfbba;
        (dist.f[dirB])[kb]     = mfbbc;
        (dist.f[dirNE])[k]     = mfaab;
        (dist.f[dirSW])[ksw]   = mfccb;
        (dist.f[dirSE])[ks]    = mfacb;
        (dist.f[dirNW])[kw]    = mfcab;
        (dist.f[dirTE])[k]     = mfaba;
        (dist.f[dirBW])[kbw]   = mfcbc;
        (dist.f[dirBE])[kb]    = mfabc;
        (dist.f[dirTW])[kw]    = mfcba;
        (dist.f[dirTN])[k]     = mfbaa;
        (dist.f[dirBS])[kbs]   = mfbcc;
        (dist.f[dirBN])[kb]    = mfbac;
        (dist.f[dirTS])[ks]    = mfbca;
        (dist.f[dirZERO])[k]   = mfbbb;
        (dist.f[dirTNE])[k]    = mfaaa;
        (dist.f[dirTSE])[ks]   = mfaca;
        (dist.f[dirBNE])[kb]   = mfaac;
        (dist.f[dirBSE])[kbs]  = mfacc;
        (dist.f[dirTNW])[kw]   = mfcaa;
        (dist.f[dirTSW])[ksw]  = mfcca;
        (dist.f[dirBNW])[kbw]  = mfcac;
        (dist.f[dirBSW])[kbsw] = mfccc;
    }
}