#ifndef LBM_CHIMERA_H
#define LBM_CHIMERA_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/Core/DataTypes.h>
#include <basics/Core/RealConstants.h>

namespace VF
{
namespace LBM
{

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardInverseChimeraWithK 
//! Transformation from distributions to central moments according to Eq. (6)-(14) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ void forwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv,
                                       real v2, real Kinverse, real K)
{
    const real m2 = mfa + mfc;
    const real m1 = mfc - mfa;
    real m0 = m2 + mfb;

    mfa = m0;
    m0 *= Kinverse;
    m0 += c1o1;
    mfb = (m1 * Kinverse - m0 * vv) * K;
    mfc = ((m2 - c2o1 * m1 * vv) * Kinverse + v2 * m0) * K;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardInverseChimeraWithK
//! Transformation from central moments to distributions according to Eq. (57)-(65) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! ] Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ void backwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv,
                                        real v2, real Kinverse, real K)
{
    const real m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 - vv) * c1o2) * K;
    const real m1 = (((mfa - mfc) - c2o1 * mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (-v2)) * K;

    mfc = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 + vv) * c1o2) * K;
    mfa = m0;
    mfb = m1;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardChimera 
//! Transformation from distributions to central moments according to Eq. (6)-(14) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations.
//! Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    const real m1 = (mfa + mfc) + mfb;
    const real m2 = mfc - mfa;

    mfc = (mfc + mfa) + (v2 * m1 - c2o1 * vv * m2);
    mfb = m2 - vv * m1;
    mfa = m1;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardChimera 
//! Transformation from central moments to distributions according to Eq. (57)-(65) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations.
//! Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
inline __host__ __device__ void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    const real ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
    const real mb = ((mfa - mfc) - mfa * v2) - c2o1 * mfb * vv;

    mfc = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
    mfb = mb;
    mfa = ma;
}


inline __host__ __device__ void forwardChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real K) 
{

    const real m2 = mfa + mfc;
    const real m1 = mfc - mfa;
    const real m0 = m2 + mfb;
    mfa = m0;
    //m0     += K;
    mfb = (m1 - K*vv) - m0 * vv;
    mfc = ((m2 - c2o1*	m1 * vv) + v2*K) + v2 * m0;
    //m0 += K;
    //mfb = m1 - m0 * vv;
    //mfc = m2 - two*	m1 * vv + v2 * m0;
}


inline __host__ __device__ void backwardChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real K) 
{
    const real  m0 = (mfc - mfb)* c1o2 + mfb * (vv)+(mfa + K) * (v2 - vv) * c1o2;
    const real m1 = (mfa - mfc) - c2o1* mfb * vv + (mfa + K) * (-v2);
    mfc = (mfc + mfb)* c1o2 + mfb * (vv)+(mfa + K) * (v2 + vv) * c1o2;
    mfa = m0;
    mfb = m1;
}

}
}
#endif
