#include "BGK.h"

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



__host__ __device__ void bgk(KernelParameter parameter)
{
    auto& distribution = parameter.distribution;
    const auto omega = parameter.omega;


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
    const real drho =
        ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
        (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
        ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb; 
    const real rho = c1o1 + drho;
    const real OOrho = c1o1 / rho;    
    const real vvx = 
        ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
        (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
        (mfcbb - mfabb)) * OOrho;
    const real vvy = 
        ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
        (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
        (mfbcb - mfbab)) * OOrho;
    const real vvz = 
        ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
        (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
        (mfbbc - mfbba)) * OOrho;


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //BGK comp
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const real cusq = c3o2*(vvx*vvx + vvy*vvy + vvz*vvz);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mfbbb = mfbbb *(c1o1 + (-omega)) - (-omega)*   c8o27*  (drho - rho * cusq);
    mfcbb = mfcbb    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(vvx)+c9o2*(vvx)*(vvx)-cusq));
    mfabb = mfabb    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(-vvx) + c9o2*(-vvx)*(-vvx) - cusq));
    mfbcb = mfbcb    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(vvy)+c9o2*(vvy)*(vvy)-cusq));
    mfbab = mfbab    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(-vvy) + c9o2*(-vvy)*(-vvy) - cusq));
    mfbbc = mfbbc    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(vvz)+c9o2*(vvz)*(vvz)-cusq));
    mfbba = mfbba    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + rho * (c3o1*(-vvz) + c9o2*(-vvz)*(-vvz) - cusq));
    mfccb = mfccb   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + rho * (c3o1*(vvx + vvy) + c9o2*(vvx + vvy)*(vvx + vvy) - cusq));
    mfaab = mfaab   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + rho * (c3o1*(-vvx - vvy) + c9o2*(-vvx - vvy)*(-vvx - vvy) - cusq));
    mfcab = mfcab   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vvx - vvy) + c9o2*(vvx - vvy)*(vvx - vvy) - cusq));
    mfacb = mfacb   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vvx + vvy) + c9o2*(-vvx + vvy)*(-vvx + vvy) - cusq));
    mfcbc = mfcbc   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vvx + vvz) + c9o2*(vvx + vvz)*(vvx + vvz) - cusq));
    mfaba = mfaba   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vvx - vvz) + c9o2*(-vvx - vvz)*(-vvx - vvz) - cusq));
    mfcba = mfcba   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vvx - vvz) + c9o2*(vvx - vvz)*(vvx - vvz) - cusq));
    mfabc = mfabc   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vvx + vvz) + c9o2*(-vvx + vvz)*(-vvx + vvz) - cusq));
    mfbcc = mfbcc   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vvy + vvz) + c9o2*(vvy + vvz)*(vvy + vvz) - cusq));
    mfbaa = mfbaa   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vvy - vvz) + c9o2*(-vvy - vvz)*(-vvy - vvz) - cusq));
    mfbca = mfbca   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(vvy - vvz) + c9o2*(vvy - vvz)*(vvy - vvz) - cusq));
    mfbac = mfbac   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + rho * (c3o1*(-vvy + vvz) + c9o2*(-vvy + vvz)*(-vvy + vvz) - cusq));
    mfccc = mfccc  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vvx + vvy + vvz) + c9o2*(vvx + vvy + vvz)*(vvx + vvy + vvz) - cusq));
    mfaaa = mfaaa  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vvx - vvy - vvz) + c9o2*(-vvx - vvy - vvz)*(-vvx - vvy - vvz) - cusq));
    mfcca = mfcca  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vvx + vvy - vvz) + c9o2*(vvx + vvy - vvz)*(vvx + vvy - vvz) - cusq));
    mfaac = mfaac  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vvx - vvy + vvz) + c9o2*(-vvx - vvy + vvz)*(-vvx - vvy + vvz) - cusq));
    mfcac = mfcac  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vvx - vvy + vvz) + c9o2*(vvx - vvy + vvz)*(vvx - vvy + vvz) - cusq));
    mfaca = mfaca  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vvx + vvy - vvz) + c9o2*(-vvx + vvy - vvz)*(-vvx + vvy - vvz) - cusq));
    mfcaa = mfcaa  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(vvx - vvy - vvz) + c9o2*(vvx - vvy - vvz)*(vvx - vvy - vvz) - cusq));
    mfacc = mfacc  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + rho * (c3o1*(-vvx + vvy + vvz) + c9o2*(-vvx + vvy + vvz)*(-vvx + vvy + vvz) - cusq));

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

