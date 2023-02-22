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
//! \file ChimeraTransformation.h
//! \ingroup LBM/GPUHelperFunctions
//! \author Martin Schoenherr, Anna Wellmann, Soeren Peters
//=======================================================================================
#ifndef CHIMERA_TRANSFORMATION_H
#define CHIMERA_TRANSFORMATION_H

#include "LBM/LB.h"

#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

namespace vf::gpu
{

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardInverseChimeraWithK
//! Transformation from distributions to central moments according to Eq. (6)-(14) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> Modified for lower round-off errors.
__inline__ __device__ void forwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K)
{
    real m2 = mfa + mfc;
    real m1 = mfc - mfa;
    real m0 = m2 + mfb;
    mfa = m0;
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
__inline__ __device__ void backwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv, real v2, real Kinverse, real K)
{
    real m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 - vv) * c1o2) * K;
    real m1 = (((mfa - mfc) - c2o1 * mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (-v2)) * K;
    mfc = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + c1o1) * (v2 + vv) * c1o2) * K;
    mfa = m0;
    mfb = m1;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardChimera
//! Transformation from distributions to central moments according to Eq. (6)-(14) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations. Modified for lower round-off
//! errors.
__inline__ __device__ void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    real m1 = (mfa + mfc) + mfb;
    real m2 = mfc - mfa;
    mfc = (mfc + mfa) + (v2 * m1 - c2o1 * vv * m2);
    mfb = m2 - vv * m1;
    mfa = m1;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardChimera
//! Transformation from central moments to distributions according to Eq. (57)-(65) in \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
//! ]</b></a> for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations. Modified for lower round-off
//! errors.
__inline__ __device__ void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2)
{
    real ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
    real mb = ((mfa - mfc) - mfa * v2) - c2o1 * mfb * vv;
    mfc = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
    mfb = mb;
    mfa = ma;
}

} // namespace vf::gpu

#endif