#include "BackwardChimera.h"

#include <basics/Core/RealConstants.h>


void VF::LBM::forwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv,
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

void VF::LBM::backwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv,
                                        real v2, real Kinverse, real K)
{
    const real m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 - vv) * c1o2) * K;
    const real m1 = (((mfa - mfc) - c2o1 * mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (-v2)) * K;

    mfc = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 + vv) * c1o2) * K;
    mfa = m0;
    mfb = m1;
}

void VF::LBM::forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    const real m1 = (mfa + mfc) + mfb;
    const real m2 = mfc - mfa;

    mfc = (mfc + mfa) + (v2 * m1 - c2o1 * vv * m2);
    mfb = m2 - vv * m1;
    mfa = m1;
}

void VF::LBM::backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    const real ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
    const real mb = ((mfa - mfc) - mfa * v2) - c2o1 * mfb * vv;

    mfc = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
    mfb = mb;
    mfa = ma;
}
