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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Data Data
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "EsoSplit.h"

EsoSplit::EsoSplit() = default;
//////////////////////////////////////////////////////////////////////////
EsoSplit::EsoSplit(size_t nx1, size_t nx2, size_t nx3, real value)
{
    this->NX1 = nx1;
    this->NX2 = nx2;
    this->NX3 = nx3;

    this->localDistributions =
        std::make_shared<CbArray4D<real, IndexerX4X3X2X1>>(13, nx1 + 1, nx2 + 1, nx3 + 1, value);
    this->nonLocalDistributions =
        std::make_shared<CbArray4D<real, IndexerX4X3X2X1>>(13, nx1 + 1, nx2 + 1, nx3 + 1, value);

    this->zeroDistributions = std::make_shared<CbArray3D<real, IndexerX3X2X1>>(nx1, nx2, nx3, value);
}
//////////////////////////////////////////////////////////////////////////
EsoSplit::~EsoSplit() = default;
//////////////////////////////////////////////////////////////////////////
void EsoSplit::swap() { std::swap(this->localDistributions, this->nonLocalDistributions); }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    const size_t x1p = x1 + 1;
    const size_t x2p = x2 + 1;
    const size_t x3p = x3 + 1;

    f[dP00] = (*this->localDistributions)(eP00, x1, x2, x3);
    f[d0P0] = (*this->localDistributions)(e0P0, x1, x2, x3);
    f[d00P] = (*this->localDistributions)(e00P, x1, x2, x3);
    f[dPP0] = (*this->localDistributions)(ePP0, x1, x2, x3);
    f[dMP0] = (*this->localDistributions)(eMP0, x1p, x2, x3);
    f[dP0P] = (*this->localDistributions)(eP0P, x1, x2, x3);
    f[dM0P] = (*this->localDistributions)(eM0P, x1p, x2, x3);
    f[d0PP] = (*this->localDistributions)(e0PP, x1, x2, x3);
    f[d0MP] = (*this->localDistributions)(e0MP, x1, x2p, x3);
    f[dPPP] = (*this->localDistributions)(ePPP, x1, x2, x3);
    f[dMPP] = (*this->localDistributions)(eMPP, x1p, x2, x3);
    f[dPMP] = (*this->localDistributions)(ePMP, x1, x2p, x3);
    f[dMMP] = (*this->localDistributions)(eMMP, x1p, x2p, x3);

    f[dM00] = (*this->nonLocalDistributions)(eM00, x1p, x2, x3);
    f[d0M0] = (*this->nonLocalDistributions)(e0M0, x1, x2p, x3);
    f[d00M] = (*this->nonLocalDistributions)(e00M, x1, x2, x3p);
    f[dMM0] = (*this->nonLocalDistributions)(eMM0, x1p, x2p, x3);
    f[dPM0] = (*this->nonLocalDistributions)(ePM0, x1, x2p, x3);
    f[dM0M] = (*this->nonLocalDistributions)(eM0M, x1p, x2, x3p);
    f[dP0M] = (*this->nonLocalDistributions)(eP0M, x1, x2, x3p);
    f[d0MM] = (*this->nonLocalDistributions)(e0MM, x1, x2p, x3p);
    f[d0PM] = (*this->nonLocalDistributions)(e0PM, x1, x2, x3p);
    f[dMMM] = (*this->nonLocalDistributions)(eMMM, x1p, x2p, x3p);
    f[dPMM] = (*this->nonLocalDistributions)(ePMM, x1, x2p, x3p);
    f[dMPM] = (*this->nonLocalDistributions)(eMPM, x1p, x2, x3p);
    f[dPPM] = (*this->nonLocalDistributions)(ePPM, x1, x2, x3p);

    f[d000] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    (*this->localDistributions)(eP00, x1, x2, x3) = f[iP00];
    (*this->localDistributions)(e0P0, x1, x2, x3) = f[i0P0];
    (*this->localDistributions)(e00P, x1, x2, x3) = f[i00P];
    (*this->localDistributions)(ePP0, x1, x2, x3) = f[iPP0];
    (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = f[iMP0];
    (*this->localDistributions)(eP0P, x1, x2, x3) = f[iP0P];
    (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = f[iM0P];
    (*this->localDistributions)(e0PP, x1, x2, x3) = f[i0PP];
    (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = f[i0MP];
    (*this->localDistributions)(ePPP, x1, x2, x3) = f[iPPP];
    (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = f[iMPP];
    (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = f[iPMP];
    (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = f[iMMP];

    (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = f[iM00];
    (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = f[i0M0];
    (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = f[i00M];
    (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = f[iMM0];
    (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = f[iPM0];
    (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = f[iM0M];
    (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = f[iP0M];
    (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = f[i0MM];
    (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = f[i0PM];
    (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[iMMM];
    (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = f[iPMM];
    (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = f[iMPM];
    (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = f[iPPM];

    (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    f[iP00] = (*this->localDistributions)(eP00, x1, x2, x3);
    f[i0P0] = (*this->localDistributions)(e0P0, x1, x2, x3);
    f[i00P] = (*this->localDistributions)(e00P, x1, x2, x3);
    f[iPP0] = (*this->localDistributions)(ePP0, x1, x2, x3);
    f[iMP0] = (*this->localDistributions)(eMP0, x1 + 1, x2, x3);
    f[iP0P] = (*this->localDistributions)(eP0P, x1, x2, x3);
    f[iM0P] = (*this->localDistributions)(eM0P, x1 + 1, x2, x3);
    f[i0PP] = (*this->localDistributions)(e0PP, x1, x2, x3);
    f[i0MP] = (*this->localDistributions)(e0MP, x1, x2 + 1, x3);
    f[iPPP] = (*this->localDistributions)(ePPP, x1, x2, x3);
    f[iMPP] = (*this->localDistributions)(eMPP, x1 + 1, x2, x3);
    f[iPMP] = (*this->localDistributions)(ePMP, x1, x2 + 1, x3);
    f[iMMP] = (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3);

    f[iM00] = (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3);
    f[i0M0] = (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3);
    f[i00M] = (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1);
    f[iMM0] = (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3);
    f[iPM0] = (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3);
    f[iM0M] = (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1);
    f[iP0M] = (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1);
    f[i0MM] = (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1);
    f[i0PM] = (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1);
    f[iMMM] = (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
    f[iPMM] = (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1);
    f[iMPM] = (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1);
    f[iPPM] = (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1);

    f[d000] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    (*this->localDistributions)(eP00, x1, x2, x3) = f[dP00];
    (*this->localDistributions)(e0P0, x1, x2, x3) = f[d0P0];
    (*this->localDistributions)(e00P, x1, x2, x3) = f[d00P];
    (*this->localDistributions)(ePP0, x1, x2, x3) = f[dPP0];
    (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = f[dMP0];
    (*this->localDistributions)(eP0P, x1, x2, x3) = f[dP0P];
    (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = f[dM0P];
    (*this->localDistributions)(e0PP, x1, x2, x3) = f[d0PP];
    (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = f[d0MP];
    (*this->localDistributions)(ePPP, x1, x2, x3) = f[dPPP];
    (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = f[dMPP];
    (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = f[dPMP];
    (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = f[dMMP];

    (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = f[dM00];
    (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = f[d0M0];
    (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = f[d00M];
    (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = f[dMM0];
    (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = f[dPM0];
    (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = f[dM0M];
    (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = f[dP0M];
    (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = f[d0MM];
    (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = f[d0PM];
    (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[dMMM];
    (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = f[dPMM];
    (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = f[dMPM];
    (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = f[dPPM];

    (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                                unsigned long int direction)
{
    using namespace vf::lbm::dir;

    if ((direction & etP00) == etP00)
        (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = f[dP00];
    if ((direction & etM00) == etM00)
        (*this->localDistributions)(eP00, x1, x2, x3) = f[dM00];
    if ((direction & et0M0) == et0M0)
        (*this->localDistributions)(e0P0, x1, x2, x3) = f[d0M0];
    if ((direction & et0P0) == et0P0)
        (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = f[d0P0];
    if ((direction & et00M) == et00M)
        (*this->localDistributions)(e00P, x1, x2, x3) = f[d00M];
    if ((direction & et00P) == et00P)
        (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = f[d00P];
    if ((direction & etMM0) == etMM0)
        (*this->localDistributions)(ePP0, x1, x2, x3) = f[dMM0];
    if ((direction & etPP0) == etPP0)
        (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = f[dPP0];
    if ((direction & etMP0) == etMP0)
        (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = f[dMP0];
    if ((direction & etPM0) == etPM0)
        (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = f[dPM0];
    if ((direction & etM0M) == etM0M)
        (*this->localDistributions)(eP0P, x1, x2, x3) = f[dM0M];
    if ((direction & etP0P) == etP0P)
        (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = f[dP0P];
    if ((direction & etM0P) == etM0P)
        (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = f[dM0P];
    if ((direction & etP0M) == etP0M)
        (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = f[dP0M];
    if ((direction & et0MM) == et0MM)
        (*this->localDistributions)(e0PP, x1, x2, x3) = f[d0MM];
    if ((direction & et0PP) == et0PP)
        (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = f[d0PP];
    if ((direction & et0MP) == et0MP)
        (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = f[d0MP];
    if ((direction & et0PM) == et0PM)
        (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = f[d0PM];
    if ((direction & etMMM) == etMMM)
        (*this->localDistributions)(ePPP, x1, x2, x3) = f[dMMM];
    if ((direction & etPPP) == etPPP)
        (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[dPPP];
    if ((direction & etPMM) == etPMM)
        (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = f[dPMM];
    if ((direction & etMPP) == etMPP)
        (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = f[dMPP];
    if ((direction & etMPM) == etMPM)
        (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = f[dMPM];
    if ((direction & etPMP) == etPMP)
        (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = f[dPMP];
    if ((direction & etPPM) == etPPM)
        (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = f[dPPM];
    if ((direction & etMMP) == etMMP)
        (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = f[dMMP];
    if ((direction & et000) == et000)
        (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                                int direction)
{
    using namespace vf::lbm::dir;
 
    switch (direction) {
        case dP00:
            (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = f;
            break;
        case dM00:
            (*this->localDistributions)(eP00, x1, x2, x3) = f;
            break;
        case d0M0:
            (*this->localDistributions)(e0P0, x1, x2, x3) = f;
            break;
        case d0P0:
            (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = f;
            break;
        case d00M:
            (*this->localDistributions)(e00P, x1, x2, x3) = f;
            break;
        case d00P:
            (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = f;
            break;
        case dMM0:
            (*this->localDistributions)(ePP0, x1, x2, x3) = f;
            break;
        case dPP0:
            (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = f;
            break;
        case dMP0:
            (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = f;
            break;
        case dPM0:
            (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = f;
            break;
        case dM0M:
            (*this->localDistributions)(eP0P, x1, x2, x3) = f;
            break;
        case dP0P:
            (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = f;
            break;
        case dM0P:
            (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = f;
            break;
        case dP0M:
            (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = f;
            break;
        case d0MM:
            (*this->localDistributions)(e0PP, x1, x2, x3) = f;
            break;
        case d0PP:
            (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = f;
            break;
        case d0MP:
            (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = f;
            break;
        case d0PM:
            (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = f;
            break;
        case dMMM:
            (*this->localDistributions)(ePPP, x1, x2, x3) = f;
            break;
        case dPPP:
            (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f;
            break;
        case dPMM:
            (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = f;
            break;
        case dMPP:
            (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = f;
            break;
        case dMPM:
            (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = f;
            break;
        case dPMP:
            (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = f;
            break;
        case dPPM:
            (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = f;
            break;
        case dMMP:
            (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = f;
            break;
        case d000:
            (*this->zeroDistributions)(x1, x2, x3) = f;
            break;
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2,
                                                                   size_t x3, unsigned long int direction)
{
    using namespace vf::lbm::dir;

    if ((direction & etP00) == etP00)
        (*this->localDistributions)(eP00, x1, x2, x3) = f[dP00];
    if ((direction & etM00) == etM00)
        (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = f[dM00];
    if ((direction & et0M0) == et0M0)
        (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = f[d0M0];
    if ((direction & et0P0) == et0P0)
        (*this->localDistributions)(e0P0, x1, x2, x3) = f[d0P0];
    if ((direction & et00M) == et00M)
        (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = f[d00M];
    if ((direction & et00P) == et00P)
        (*this->localDistributions)(e00P, x1, x2, x3) = f[d00P];
    if ((direction & etMM0) == etMM0)
        (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = f[dMM0];
    if ((direction & etPP0) == etPP0)
        (*this->localDistributions)(ePP0, x1, x2, x3) = f[dPP0];
    if ((direction & etMP0) == etMP0)
        (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = f[dMP0];
    if ((direction & etPM0) == etPM0)
        (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = f[dPM0];
    if ((direction & etM0M) == etM0M)
        (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = f[dM0M];
    if ((direction & etP0P) == etP0P)
        (*this->localDistributions)(eP0P, x1, x2, x3) = f[dP0P];
    if ((direction & etM0P) == etM0P)
        (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = f[dM0P];
    if ((direction & etP0M) == etP0M)
        (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = f[dP0M];
    if ((direction & et0MM) == et0MM)
        (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = f[d0MM];
    if ((direction & et0PP) == et0PP)
        (*this->localDistributions)(e0PP, x1, x2, x3) = f[d0PP];
    if ((direction & et0MP) == et0MP)
        (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = f[d0MP];
    if ((direction & et0PM) == et0PM)
        (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = f[d0PM];
    if ((direction & etMMM) == etMMM)
        (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[dMMM];
    if ((direction & etPPP) == etPPP)
        (*this->localDistributions)(ePPP, x1, x2, x3) = f[dPPP];
    if ((direction & etPMM) == etPMM)
        (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = f[dPMM];
    if ((direction & etMPP) == etMPP)
        (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = f[dMPP];
    if ((direction & etMPM) == etMPM)
        (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = f[dMPM];
    if ((direction & etPMP) == etPMP)
        (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = f[dPMP];
    if ((direction & etPPM) == etPPM)
        (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = f[dPPM];
    if ((direction & etMMP) == etMMP)
        (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = f[dMMP];
    if ((direction & et000) == et000)
        (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                                   unsigned long int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dP00:
            (*this->localDistributions)(eP00, x1, x2, x3) = f;
            break;
        case dM00:
            (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = f;
            break;
        case d0M0:
            (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = f;
            break;
        case d0P0:
            (*this->localDistributions)(e0P0, x1, x2, x3) = f;
            break;
        case d00M:
            (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = f;
            break;
        case d00P:
            (*this->localDistributions)(e00P, x1, x2, x3) = f;
            break;
        case dMM0:
            (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = f;
            break;
        case dPP0:
            (*this->localDistributions)(ePP0, x1, x2, x3) = f;
            break;
        case dMP0:
            (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = f;
            break;
        case dPM0:
            (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = f;
            break;
        case dM0M:
            (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = f;
            break;
        case dP0P:
            (*this->localDistributions)(eP0P, x1, x2, x3) = f;
            break;
        case dM0P:
            (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = f;
            break;
        case dP0M:
            (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = f;
            break;
        case d0MM:
            (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = f;
            break;
        case d0PP:
            (*this->localDistributions)(e0PP, x1, x2, x3) = f;
            break;
        case d0MP:
            (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = f;
            break;
        case d0PM:
            (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = f;
            break;
        case dMMM:
            (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f;
            break;
        case dPPP:
            (*this->localDistributions)(ePPP, x1, x2, x3) = f;
            break;
        case dPMM:
            (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = f;
            break;
        case dMPP:
            (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = f;
            break;
        case dMPM:
            (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = f;
            break;
        case dPMP:
            (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = f;
            break;
        case dPPM:
            (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = f;
            break;
        case dMMP:
            (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = f;
            break;
        case d000:
            (*this->zeroDistributions)(x1, x2, x3) = f;
            break;
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
real EsoSplit::getPreCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dM00:
            return (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3);
        case dP00:
            return (*this->localDistributions)(eP00, x1, x2, x3);
        case d0P0:
            return (*this->localDistributions)(e0P0, x1, x2, x3);
        case d0M0:
            return (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3);
        case d00P:
            return (*this->localDistributions)(e00P, x1, x2, x3);
        case d00M:
            return (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1);
        case dPP0:
            return (*this->localDistributions)(ePP0, x1, x2, x3);
        case dMM0:
            return (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3);
        case dPM0:
            return (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3);
        case dMP0:
            return (*this->localDistributions)(eMP0, x1 + 1, x2, x3);
        case dP0P:
            return (*this->localDistributions)(eP0P, x1, x2, x3);
        case dM0M:
            return (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1);
        case dP0M:
            return (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1);
        case dM0P:
            return (*this->localDistributions)(eM0P, x1 + 1, x2, x3);
        case d0PP:
            return (*this->localDistributions)(e0PP, x1, x2, x3);
        case d0MM:
            return (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1);
        case d0PM:
            return (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1);
        case d0MP:
            return (*this->localDistributions)(e0MP, x1, x2 + 1, x3);
        case dPPP:
            return (*this->localDistributions)(ePPP, x1, x2, x3);
        case dMMM:
            return (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
        case dMPP:
            return (*this->localDistributions)(eMPP, x1 + 1, x2, x3);
        case dPMM:
            return (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1);
        case dPMP:
            return (*this->localDistributions)(ePMP, x1, x2 + 1, x3);
        case dMPM:
            return (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1);
        case dMMP:
            return (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3);
        case dPPM:
            return (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1);
        case d000:
            return (*this->zeroDistributions)(x1, x2, x3);
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
real EsoSplit::getPostCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dP00:
            return (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3);
        case dM00:
            return (*this->localDistributions)(eP00, x1, x2, x3);
        case d0M0:
            return (*this->localDistributions)(e0P0, x1, x2, x3);
        case d0P0:
            return (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3);
        case d00M:
            return (*this->localDistributions)(e00P, x1, x2, x3);
        case d00P:
            return (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1);
        case dMM0:
            return (*this->localDistributions)(ePP0, x1, x2, x3);
        case dPP0:
            return (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3);
        case dMP0:
            return (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3);
        case dPM0:
            return (*this->localDistributions)(eMP0, x1 + 1, x2, x3);
        case dM0M:
            return (*this->localDistributions)(eP0P, x1, x2, x3);
        case dP0P:
            return (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1);
        case dM0P:
            return (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1);
        case dP0M:
            return (*this->localDistributions)(eM0P, x1 + 1, x2, x3);
        case d0MM:
            return (*this->localDistributions)(e0PP, x1, x2, x3);
        case d0PP:
            return (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1);
        case d0MP:
            return (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1);
        case d0PM:
            return (*this->localDistributions)(e0MP, x1, x2 + 1, x3);
        case dMMM:
            return (*this->localDistributions)(ePPP, x1, x2, x3);
        case dPPP:
            return (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
        case dPMM:
            return (*this->localDistributions)(eMPP, x1 + 1, x2, x3);
        case dMPP:
            return (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1);
        case dMPM:
            return (*this->localDistributions)(ePMP, x1, x2 + 1, x3);
        case dPMP:
            return (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1);
        case dPPM:
            return (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3);
        case dMMP:
            return (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1);
        case d000:
            return (*this->zeroDistributions)(x1, x2, x3);
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
size_t EsoSplit::getNX1() const { return NX1; }
//////////////////////////////////////////////////////////////////////////
size_t EsoSplit::getNX2() const { return NX2; }
//////////////////////////////////////////////////////////////////////////
size_t EsoSplit::getNX3() const { return NX3; }
//////////////////////////////////////////////////////////////////////////
CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr EsoSplit::getLocalDistributions()
{
    return this->localDistributions;
}
//////////////////////////////////////////////////////////////////////////
CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr EsoSplit::getNonLocalDistributions()
{
    return this->nonLocalDistributions;
}
//////////////////////////////////////////////////////////////////////////
CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr EsoSplit::getZeroDistributions()
{
    return this->zeroDistributions;
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNX1(size_t newNX1) { NX1 = newNX1; }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNX2(size_t newNX2) { NX2 = newNX2; }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNX3(size_t newNX3) { NX3 = newNX3; }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array)
{
    localDistributions = array;
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNonLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array)
{
    nonLocalDistributions = array;
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setZeroDistributions(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr array)
{
    zeroDistributions = array;
}

//////////////////////////////////////////////////////////////////////////

//! \}
