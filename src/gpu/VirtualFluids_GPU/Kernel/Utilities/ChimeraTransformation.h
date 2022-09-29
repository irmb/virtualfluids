#ifndef CHIMERA_TRANSFORMATION_H
#define CHIMERA_TRANSFORMATION_H

#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardInverseChimeraWithK
//! Transformation from distributions to central moments according to Eq. (6)-(14) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> Modified for lower round-off errors.
inline __device__ void forwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K)
{
    real m2 = mfa + mfc;
    real m1 = mfc - mfa;
    real m0 = m2 + mfb;
    mfa     = m0;
    m0 *= Kinverse;
    m0 += c1o1;
    mfb = (m1 * Kinverse - m0 * vv) * K;
    mfc = ((m2 - c2o1 * m1 * vv) * Kinverse + v2 * m0) * K;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardInverseChimeraWithK
//! Transformation from central moments to distributions according to Eq. (57)-(65) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> Modified for lower round-off errors.
inline __device__ void backwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K)
{
    real m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 - vv) * c1o2) * K;
    real m1 = (((mfa - mfc) - c2o1 * mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (-v2)) * K;
    mfc     = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 + vv) * c1o2) * K;
    mfa     = m0;
    mfb     = m1;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardChimera
//! Transformation from distributions to central moments according to Eq. (6)-(14) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations. Modified for lower round-off
//! errors.
inline __device__ void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    real m1 = (mfa + mfc) + mfb;
    real m2 = mfc - mfa;
    mfc     = (mfc + mfa) + (v2 * m1 - c2o1 * vv * m2);
    mfb     = m2 - vv * m1;
    mfa     = m1;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardChimera
//! Transformation from central moments to distributions according to Eq. (57)-(65) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations. Modified for lower round-off
//! errors.
inline __device__ void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    real ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
    real mb = ((mfa - mfc) - mfa * v2) - c2o1 * mfb * vv;
    mfc     = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
    mfb     = mb;
    mfa     = ma;
}
#endif