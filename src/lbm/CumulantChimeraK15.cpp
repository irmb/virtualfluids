#include "CumulantChimeraK17.h"

#include <cmath>

#include <basics/Core/DataTypes.h>
#include <basics/Core/RealConstants.h>

#include "constants/NumericConstants.h"
#include "constants/D3Q27.h"

#include "Chimera.h"
#include "MacroscopicQuantities.h"

namespace vf
{
namespace lbm
{

using namespace constant;


//////////////////////////////////////////////////////////////////////////
//! Cumulant K17 Kernel is based on \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! and \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
//////////////////////////////////////////////////////////////////////////
__host__ __device__ void cumulantChimeraK15(Distribution27& distribution, real omega, real* forces)
{
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to 
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
    //!
    real mfcbb = distribution.f[dir::PZZ];
    real mfabb = distribution.f[dir::MZZ];
    real mfbcb = distribution.f[dir::ZPZ];
    real mfbab = distribution.f[dir::ZMZ];
    real mfbbc = distribution.f[dir::ZZP];
    real mfbba = distribution.f[dir::ZZM];
    real mfccb = distribution.f[dir::PPZ];
    real mfaab = distribution.f[dir::MMZ];
    real mfcab = distribution.f[dir::PMZ];
    real mfacb = distribution.f[dir::MPZ];
    real mfcbc = distribution.f[dir::PZP];
    real mfaba = distribution.f[dir::MZM];
    real mfcba = distribution.f[dir::PZM];
    real mfabc = distribution.f[dir::MZP];
    real mfbcc = distribution.f[dir::ZPP];
    real mfbaa = distribution.f[dir::ZMM];
    real mfbca = distribution.f[dir::ZPM];
    real mfbac = distribution.f[dir::ZMP];
    real mfccc = distribution.f[dir::PPP];
    real mfacc = distribution.f[dir::MPP];
    real mfcac = distribution.f[dir::PMP];
    real mfaac = distribution.f[dir::MMP];
    real mfcca = distribution.f[dir::PPM];
    real mfaca = distribution.f[dir::MPM];
    real mfcaa = distribution.f[dir::PMM];
    real mfaaa = distribution.f[dir::MMM];
    real mfbbb = distribution.f[dir::ZZZ];

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
    //!
    real drho =
        ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
        (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
        ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb; 
    real rho = c1o1 + drho;
    real OOrho = c1o1 / rho;    
    real vvx = 
        ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
        (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
        (mfcbb - mfabb)) * OOrho;
    real vvy = 
        ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
        (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
        (mfbcb - mfbab)) * OOrho;
    real vvz = 
        ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
        (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
        (mfbbc - mfbba)) * OOrho;
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Add half of the acceleration (body force) to the velocity as in Eq. (42) \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
    //!
    vvx += forces[0] * c1o2;
    vvy += forces[1] * c1o2;
    vvz += forces[2] * c1o2;
    ////////////////////////////////////////////////////////////////////////////////////
    // calculate the square of velocities for this lattice node
    real vx2 = vvx*vvx;
    real vy2 = vvy*vvy;
    real vz2 = vvz*vvz;
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set relaxation limiters for third order cumulants to default value \f$ \lambda=0.001 \f$ according to section 6 in \ref
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    real wadjust;
    real qudricLimitP = c1o100;
    real qudricLimitM = c1o100;
    real qudricLimitD = c1o100;
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from well conditioned distributions to central moments as defined in Appendix J in \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
    //! see also Eq. (6)-(14) in \ref
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    vf::lbm::forwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36);
    vf::lbm::forwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::forwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36);
    vf::lbm::forwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::forwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2,  c9o4,  c4o9);
    vf::lbm::forwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::forwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36);
    vf::lbm::forwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::forwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36);   
    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    vf::lbm::forwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2,  c6o1,  c1o6);
    vf::lbm::forwardChimera(            mfaab, mfabb, mfacb, vvy, vy2);
    vf::lbm::forwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18);
    vf::lbm::forwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2,  c3o2,  c2o3);
    vf::lbm::forwardChimera(            mfbab, mfbbb, mfbcb, vvy, vy2);
    vf::lbm::forwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2,  c9o2,  c2o9);
    vf::lbm::forwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2,  c6o1,  c1o6);
    vf::lbm::forwardChimera(            mfcab, mfcbb, mfccb, vvy, vy2);
    vf::lbm::forwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18);   
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    vf::lbm::forwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
    vf::lbm::forwardChimera(            mfaba, mfbba, mfcba, vvx, vx2);
    vf::lbm::forwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3);
    vf::lbm::forwardChimera(            mfaab, mfbab, mfcab, vvx, vx2);
    vf::lbm::forwardChimera(            mfabb, mfbbb, mfcbb, vvx, vx2);
    vf::lbm::forwardChimera(            mfacb, mfbcb, mfccb, vvx, vx2);
    vf::lbm::forwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3);
    vf::lbm::forwardChimera(            mfabc, mfbbc, mfcbc, vvx, vx2);
    vf::lbm::forwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c3o1, c1o9); 
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Setting relaxation rates for non-hydrodynamic cumulants (default values). Variable names and equations    according to
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!  => [NAME IN PAPER]=[NAME IN CODE]=[DEFAULT VALUE].
    //!  - Trace of second order cumulants \f$ C_{200}+C_{020}+C_{002} \f$ used to adjust bulk  viscosity:\f$\omega_2=OxxPyyPzz=1.0 \f$.
    //!  - Third order cumulants \f$ C_{120}+C_{102}, C_{210}+C_{012}, C_{201}+C_{021} \f$: \f$ \omega_3=OxyyPxzz   \f$ set according to Eq. (111) with simplifications assuming \f$ \omega_2=1.0\f$.
    //!  - Third order cumulants \f$ C_{120}-C_{102}, C_{210}-C_{012}, C_{201}-C_{021} \f$: \f$ \omega_4 =  OxyyMxzz \f$ set according to Eq. (112) with simplifications assuming \f$ \omega_2 = 1.0\f$.
    //!  - Third order cumulants \f$ C_{111} \f$: \f$ \omega_5 = Oxyz \f$ set according to Eq. (113) with   simplifications assuming \f$ \omega_2 = 1.0\f$  (modify for different bulk viscosity).
    //!  - Fourth order cumulants \f$ C_{220}, C_{202}, C_{022}, C_{211}, C_{121}, C_{112} \f$: for simplification  all set to the same default value \f$ \omega_6=\omega_7=\omega_8=O4=1.0 \f$.
    //!  - Fifth order cumulants \f$ C_{221}, C_{212}, C_{122}\f$: \f$\omega_9=O5=1.0\f$.
    //!  - Sixth order cumulant \f$ C_{222}\f$: \f$\omega_{10}=O6=1.0\f$.
    //!
    ////////////////////////////////////////////////////////////
    //2.
    real OxxPyyPzz = c1o1;
    ////////////////////////////////////////////////////////////
    //3.
    real OxyyPxzz = c1o1;//three  * (two - omega) / (three  - omega);//one;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
    real OxyyMxzz = c1o1;//six    * (two - omega) / (six    - omega);//one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
    real Oxyz = c1o1;//twelve * (two - omega) / (twelve + omega);//one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
    ////////////////////////////////////////////////////////////
    //4.
    real O4 = c1o1;
    ////////////////////////////////////////////////////////////
    //5.
    real O5 = c1o1;
    ////////////////////////////////////////////////////////////
    //6.
    real O6 = c1o1; 
    ////////////////////////////////////////////////////////////////////////////////////
    //! - A and B: parameters for fourth order convergence of the diffusion term according to Eq. (114) and (115) 
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //! with simplifications assuming \f$ \omega_2 = 1.0 \f$ (modify for different bulk viscosity).
    //!
    real A = (c4o1 + c2o1*omega - c3o1*omega*omega) / (c2o1 - c7o1*omega + c5o1*omega*omega);
    real B = (c4o1 + c28o1*omega - c14o1*omega*omega) / (c6o1 - c21o1*omega + c15o1*omega*omega);   
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute cumulants from central moments according to Eq. (20)-(23) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    ////////////////////////////////////////////////////////////
    //4.
    real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) * OOrho;
    real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) * OOrho;
    real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) * OOrho;  
    real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) * OOrho - c1o9*(drho   * OOrho));
    real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) * OOrho - c1o9*(drho   * OOrho));
    real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) * OOrho - c1o9*(drho   * OOrho));
    ////////////////////////////////////////////////////////////
    //5.
    real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba *  mfabc)) + c1o3 * (mfbca + mfbac)) * OOrho;
    real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba *  mfbac)) + c1o3 * (mfcba + mfabc)) * OOrho;
    real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb *  mfcba)) + c1o3 * (mfacb + mfcab)) * OOrho;
    ////////////////////////////////////////////////////////////
    //6.
    real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
        - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
        - c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
        - c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) * OOrho
        + (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
        + c2o1 * (mfcaa * mfaca * mfaac)
        + c16o1 *  mfbba * mfbab * mfabb) * OOrho * OOrho
        - c1o3 * (mfacc + mfcac + mfcca) * OOrho
        - c1o9 * (mfcaa + mfaca + mfaac) * OOrho
        + (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
        + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) * OOrho * OOrho  * c2o3
        + c1o27*((drho * drho - drho) * OOrho * OOrho));    
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute linear combinations of second and third order cumulants
    //!
    ////////////////////////////////////////////////////////////
    //2.
    real mxxPyyPzz = mfcaa + mfaca + mfaac;
    real mxxMyy = mfcaa - mfaca;
    real mxxMzz = mfcaa - mfaac;
    ////////////////////////////////////////////////////////////
    //3.
    real mxxyPyzz = mfcba + mfabc;
    real mxxyMyzz = mfcba - mfabc;  
    real mxxzPyyz = mfcab + mfacb;
    real mxxzMyyz = mfcab - mfacb;  
    real mxyyPxzz = mfbca + mfbac;
    real mxyyMxzz = mfbca - mfbac;  
    ////////////////////////////////////////////////////////////////////////////////////
    //incl. correction
    ////////////////////////////////////////////////////////////
    //! - Compute velocity  gradients from second order cumulants according to Eq. (27)-(32)
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //! Further explanations of the correction in viscosity in Appendix H of
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
    //! Note that the division by rho is omitted here as we need rho times the gradients later.
    //!
    real Dxy = -c3o1*omega*mfbba;
    real Dxz = -c3o1*omega*mfbab;
    real Dyz = -c3o1*omega*mfabb;
    real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
    real dyuy = dxux + omega * c3o2 * mxxMyy;
    real dzuz = dxux + omega * c3o2 * mxxMzz;
    ////////////////////////////////////////////////////////////
    //! - Relaxation of second order cumulants with correction terms according to Eq. (33)-(35) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2  * dzuz);
    mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
    mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);   
    ////////////////////////////////////////////////////////////////////////////////////
    ////no correction
    //mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);
    //mxxMyy += -(-omega) * (-mxxMyy);
    //mxxMzz += -(-omega) * (-mxxMzz);
    //////////////////////////////////////////////////////////////////////////
    mfabb += omega * (-mfabb);
    mfbab += omega * (-mfbab);
    mfbba += omega * (-mfbba);  
    ////////////////////////////////////////////////////////////////////////////////////
    //relax
    //////////////////////////////////////////////////////////////////////////
    // incl. limiter
    //! - Relaxation of third order cumulants including limiter according to Eq. (116)-(123)
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    wadjust   = Oxyz + (c1o1 - Oxyz)*abs_internal(mfbbb) / (abs_internal(mfbbb) + qudricLimitD);
    mfbbb    += wadjust * (-mfbbb);
    wadjust   = OxyyPxzz + (c1o1 - OxyyPxzz)*abs_internal(mxxyPyzz) / (abs_internal(mxxyPyzz) + qudricLimitP);
    mxxyPyzz += wadjust * (-mxxyPyzz);
    wadjust   = OxyyMxzz + (c1o1 - OxyyMxzz)*abs_internal(mxxyMyzz) / (abs_internal(mxxyMyzz) + qudricLimitM);
    mxxyMyzz += wadjust * (-mxxyMyzz);
    wadjust   = OxyyPxzz + (c1o1 - OxyyPxzz)*abs_internal(mxxzPyyz) / (abs_internal(mxxzPyyz) + qudricLimitP);
    mxxzPyyz += wadjust * (-mxxzPyyz);
    wadjust   = OxyyMxzz + (c1o1 - OxyyMxzz)*abs_internal(mxxzMyyz) / (abs_internal(mxxzMyyz) + qudricLimitM);
    mxxzMyyz += wadjust * (-mxxzMyyz);
    wadjust   = OxyyPxzz + (c1o1 - OxyyPxzz)*abs_internal(mxyyPxzz) / (abs_internal(mxyyPxzz) + qudricLimitP);
    mxyyPxzz += wadjust * (-mxyyPxzz);
    wadjust   = OxyyMxzz + (c1o1 - OxyyMxzz)*abs_internal(mxyyMxzz) / (abs_internal(mxyyMxzz) + qudricLimitM);
    mxyyMxzz += wadjust * (-mxyyMxzz);
    //////////////////////////////////////////////////////////////////////////
    // no limiter
    //mfbbb += OxyyMxzz * (-mfbbb);
    //mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
    //mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
    //mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
    //mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
    //mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
    //mxyyMxzz += OxyyMxzz * (-mxyyMxzz);   
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute inverse linear combinations of second and third order cumulants
    //!
    mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
    mfaca = c1o3 * (-c2o1*  mxxMyy + mxxMzz + mxxPyyPzz);
    mfaac = c1o3 * (mxxMyy - c2o1* mxxMzz + mxxPyyPzz); 
    mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
    mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
    mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
    mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
    mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
    mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
    //////////////////////////////////////////////////////////////////////////  
    //////////////////////////////////////////////////////////////////////////
    //4.
    // no limiter
    //! - Relax fourth order cumulants to modified equilibrium for fourth order convergence of diffusion according  to Eq. (43)-(48)
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    CUMacc = -O4*(c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMacc);
    CUMcac = -O4*(c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMcac);
    CUMcca = -O4*(c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (c1o1 - O4) * (CUMcca);
    CUMbbc = -O4*(c1o1 / omega - c1o2) * Dxy           * c1o3 * B + (c1o1 - O4) * (CUMbbc);
    CUMbcb = -O4*(c1o1 / omega - c1o2) * Dxz           * c1o3 * B + (c1o1 - O4) * (CUMbcb);
    CUMcbb = -O4*(c1o1 / omega - c1o2) * Dyz           * c1o3 * B + (c1o1 - O4) * (CUMcbb); 
    //////////////////////////////////////////////////////////////////////////
    //5.
    CUMbcc += O5 * (-CUMbcc);
    CUMcbc += O5 * (-CUMcbc);
    CUMccb += O5 * (-CUMccb);   
    //////////////////////////////////////////////////////////////////////////
    //6.
    CUMccc += O6 * (-CUMccc);   
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Compute central moments from post collision cumulants according to Eq. (53)-(56) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //! 
    //////////////////////////////////////////////////////////////////////////
    //4.
    mfcbb = CUMcbb + c1o3*((c3o1*mfcaa + c1o1) * mfabb + c6o1 * mfbba * mfbab) * OOrho;
    mfbcb = CUMbcb + c1o3*((c3o1*mfaca + c1o1) * mfbab + c6o1 * mfbba * mfabb) * OOrho;
    mfbbc = CUMbbc + c1o3*((c3o1*mfaac + c1o1) * mfbba + c6o1 * mfbab * mfabb) * OOrho; 
    mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba)*c9o1 + c3o1 * (mfcaa + mfaca)) * OOrho - (drho *  OOrho))*c1o9;
    mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab)*c9o1 + c3o1 * (mfcaa + mfaac)) * OOrho - (drho *  OOrho))*c1o9;
    mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb)*c9o1 + c3o1 * (mfaac + mfaca)) * OOrho - (drho *  OOrho))*c1o9; 
    //////////////////////////////////////////////////////////////////////////
    //5.
    mfbcc = CUMbcc + c1o3 *(c3o1*(mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb +    mfbba * mfabc)) + (mfbca + mfbac)) * OOrho;
    mfcbc = CUMcbc + c1o3 *(c3o1*(mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab +    mfbba * mfbac)) + (mfcba + mfabc)) * OOrho;
    mfccb = CUMccb + c1o3 *(c3o1*(mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca +    mfabb * mfcba)) + (mfacb + mfcab)) * OOrho; 
    //////////////////////////////////////////////////////////////////////////
    //6.
    mfccc =	CUMccc - ((-c4o1 *  mfbbb * mfbbb
            - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
            - c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
            - c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) * OOrho
            + (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                + c2o1 * (mfcaa * mfaca * mfaac)
                + c16o1 *  mfbba * mfbab * mfabb) * OOrho * OOrho
            - c1o3 * (mfacc + mfcac + mfcca) * OOrho
            - c1o9 * (mfcaa + mfaca + mfaac) * OOrho
            + (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) * OOrho * OOrho * c2o3
            + c1o27*((drho * drho - drho) * OOrho * OOrho));    
    ////////////////////////////////////////////////////////////////////////////////////
    //! -  Add acceleration (body force) to first order cumulants according to Eq. (85)-(87) in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
    //!
    mfbaa = -mfbaa;
    mfaba = -mfaba;
    mfaab = -mfaab; 
    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
    //! see also Eq. (88)-(96) in
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05  040 ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    vf::lbm::backwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
    vf::lbm::backwardChimera(            mfaba, mfbba, mfcba, vvx, vx2);
    vf::lbm::backwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3);
    vf::lbm::backwardChimera(            mfaab, mfbab, mfcab, vvx, vx2);
    vf::lbm::backwardChimera(            mfabb, mfbbb, mfcbb, vvx, vx2);
    vf::lbm::backwardChimera(            mfacb, mfbcb, mfccb, vvx, vx2);
    vf::lbm::backwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3);
    vf::lbm::backwardChimera(            mfabc, mfbbc, mfcbc, vvx, vx2);
    vf::lbm::backwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c9o1, c1o9);    
    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    vf::lbm::backwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2,  c6o1,  c1o6);
    vf::lbm::backwardChimera(            mfaab, mfabb, mfacb, vvy, vy2);
    vf::lbm::backwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18);
    vf::lbm::backwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2,  c3o2,  c2o3);
    vf::lbm::backwardChimera(            mfbab, mfbbb, mfbcb, vvy, vy2);
    vf::lbm::backwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2,  c9o2,  c2o9);
    vf::lbm::backwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2,  c6o1,  c1o6);
    vf::lbm::backwardChimera(            mfcab, mfcbb, mfccb, vvy, vy2);
    vf::lbm::backwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18);  
    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    vf::lbm::backwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36);
    vf::lbm::backwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::backwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36);
    vf::lbm::backwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::backwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2,  c9o4,  c4o9);
    vf::lbm::backwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::backwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36);
    vf::lbm::backwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2,  c9o1,  c1o9);
    vf::lbm::backwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36);


    ////////////////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to 
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
    //!
    distribution.f[vf::lbm::dir::MZZ] = mfcbb;
    distribution.f[vf::lbm::dir::PZZ] = mfabb;
    distribution.f[vf::lbm::dir::ZMZ] = mfbcb;
    distribution.f[vf::lbm::dir::ZPZ] = mfbab;
    distribution.f[vf::lbm::dir::ZZM] = mfbbc;
    distribution.f[vf::lbm::dir::ZZP] = mfbba;
    distribution.f[vf::lbm::dir::MMZ] = mfccb;
    distribution.f[vf::lbm::dir::PPZ] = mfaab;
    distribution.f[vf::lbm::dir::MPZ] = mfcab;
    distribution.f[vf::lbm::dir::PMZ] = mfacb;
    distribution.f[vf::lbm::dir::MZM] = mfcbc;
    distribution.f[vf::lbm::dir::PZP] = mfaba;
    distribution.f[vf::lbm::dir::MZP] = mfcba;
    distribution.f[vf::lbm::dir::PZM] = mfabc;
    distribution.f[vf::lbm::dir::ZMM] = mfbcc;
    distribution.f[vf::lbm::dir::ZPP] = mfbaa;
    distribution.f[vf::lbm::dir::ZMP] = mfbca;
    distribution.f[vf::lbm::dir::ZPM] = mfbac;
    distribution.f[vf::lbm::dir::MMM] = mfccc;
    distribution.f[vf::lbm::dir::PMM] = mfacc;
    distribution.f[vf::lbm::dir::MPM] = mfcac;
    distribution.f[vf::lbm::dir::PPM] = mfaac;
    distribution.f[vf::lbm::dir::MMP] = mfcca;
    distribution.f[vf::lbm::dir::PMP] = mfaca;
    distribution.f[vf::lbm::dir::MPP] = mfcaa;
    distribution.f[vf::lbm::dir::PPP] = mfaaa;
    distribution.f[vf::lbm::dir::ZZZ] = mfbbb;
}


}
}

