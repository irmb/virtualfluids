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

    this->splitA =
        std::make_shared<CbArray4D<real, IndexerX4X3X2X1>>(13, nx1 + 1, nx2 + 1, nx3 + 1, value);
    this->splitB =
        std::make_shared<CbArray4D<real, IndexerX4X3X2X1>>(13, nx1 + 1, nx2 + 1, nx3 + 1, value);

    this->split0 = std::make_shared<CbArray3D<real, IndexerX3X2X1>>(nx1, nx2, nx3, value);
}
//////////////////////////////////////////////////////////////////////////
EsoSplit::~EsoSplit() = default;
//////////////////////////////////////////////////////////////////////////
void EsoSplit::swap() { std::swap(this->splitA, this->splitB); }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    const size_t x1p = x1 + 1;
    const size_t x2p = x2 + 1;
    const size_t x3p = x3 + 1;

    f[dP00] = (*this->splitA)(eP00, x1, x2, x3);
    f[d0P0] = (*this->splitA)(e0P0, x1, x2, x3);
    f[d00P] = (*this->splitA)(e00P, x1, x2, x3);
    f[dPP0] = (*this->splitA)(ePP0, x1, x2, x3);
    f[dMP0] = (*this->splitA)(eMP0, x1p, x2, x3);
    f[dP0P] = (*this->splitA)(eP0P, x1, x2, x3);
    f[dM0P] = (*this->splitA)(eM0P, x1p, x2, x3);
    f[d0PP] = (*this->splitA)(e0PP, x1, x2, x3);
    f[d0MP] = (*this->splitA)(e0MP, x1, x2p, x3);
    f[dPPP] = (*this->splitA)(ePPP, x1, x2, x3);
    f[dMPP] = (*this->splitA)(eMPP, x1p, x2, x3);
    f[dPMP] = (*this->splitA)(ePMP, x1, x2p, x3);
    f[dMMP] = (*this->splitA)(eMMP, x1p, x2p, x3);

    f[dM00] = (*this->splitB)(eM00, x1p, x2, x3);
    f[d0M0] = (*this->splitB)(e0M0, x1, x2p, x3);
    f[d00M] = (*this->splitB)(e00M, x1, x2, x3p);
    f[dMM0] = (*this->splitB)(eMM0, x1p, x2p, x3);
    f[dPM0] = (*this->splitB)(ePM0, x1, x2p, x3);
    f[dM0M] = (*this->splitB)(eM0M, x1p, x2, x3p);
    f[dP0M] = (*this->splitB)(eP0M, x1, x2, x3p);
    f[d0MM] = (*this->splitB)(e0MM, x1, x2p, x3p);
    f[d0PM] = (*this->splitB)(e0PM, x1, x2, x3p);
    f[dMMM] = (*this->splitB)(eMMM, x1p, x2p, x3p);
    f[dPMM] = (*this->splitB)(ePMM, x1, x2p, x3p);
    f[dMPM] = (*this->splitB)(eMPM, x1p, x2, x3p);
    f[dPPM] = (*this->splitB)(ePPM, x1, x2, x3p);

    f[d000] = (*this->split0)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    (*this->splitA)(eP00, x1, x2, x3) = f[iP00];
    (*this->splitA)(e0P0, x1, x2, x3) = f[i0P0];
    (*this->splitA)(e00P, x1, x2, x3) = f[i00P];
    (*this->splitA)(ePP0, x1, x2, x3) = f[iPP0];
    (*this->splitA)(eMP0, x1 + 1, x2, x3) = f[iMP0];
    (*this->splitA)(eP0P, x1, x2, x3) = f[iP0P];
    (*this->splitA)(eM0P, x1 + 1, x2, x3) = f[iM0P];
    (*this->splitA)(e0PP, x1, x2, x3) = f[i0PP];
    (*this->splitA)(e0MP, x1, x2 + 1, x3) = f[i0MP];
    (*this->splitA)(ePPP, x1, x2, x3) = f[iPPP];
    (*this->splitA)(eMPP, x1 + 1, x2, x3) = f[iMPP];
    (*this->splitA)(ePMP, x1, x2 + 1, x3) = f[iPMP];
    (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = f[iMMP];

    (*this->splitB)(eM00, x1 + 1, x2, x3) = f[iM00];
    (*this->splitB)(e0M0, x1, x2 + 1, x3) = f[i0M0];
    (*this->splitB)(e00M, x1, x2, x3 + 1) = f[i00M];
    (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = f[iMM0];
    (*this->splitB)(ePM0, x1, x2 + 1, x3) = f[iPM0];
    (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = f[iM0M];
    (*this->splitB)(eP0M, x1, x2, x3 + 1) = f[iP0M];
    (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = f[i0MM];
    (*this->splitB)(e0PM, x1, x2, x3 + 1) = f[i0PM];
    (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[iMMM];
    (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = f[iPMM];
    (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = f[iMPM];
    (*this->splitB)(ePPM, x1, x2, x3 + 1) = f[iPPM];

    (*this->split0)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    f[iP00] = (*this->splitA)(eP00, x1, x2, x3);
    f[i0P0] = (*this->splitA)(e0P0, x1, x2, x3);
    f[i00P] = (*this->splitA)(e00P, x1, x2, x3);
    f[iPP0] = (*this->splitA)(ePP0, x1, x2, x3);
    f[iMP0] = (*this->splitA)(eMP0, x1 + 1, x2, x3);
    f[iP0P] = (*this->splitA)(eP0P, x1, x2, x3);
    f[iM0P] = (*this->splitA)(eM0P, x1 + 1, x2, x3);
    f[i0PP] = (*this->splitA)(e0PP, x1, x2, x3);
    f[i0MP] = (*this->splitA)(e0MP, x1, x2 + 1, x3);
    f[iPPP] = (*this->splitA)(ePPP, x1, x2, x3);
    f[iMPP] = (*this->splitA)(eMPP, x1 + 1, x2, x3);
    f[iPMP] = (*this->splitA)(ePMP, x1, x2 + 1, x3);
    f[iMMP] = (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3);

    f[iM00] = (*this->splitB)(eM00, x1 + 1, x2, x3);
    f[i0M0] = (*this->splitB)(e0M0, x1, x2 + 1, x3);
    f[i00M] = (*this->splitB)(e00M, x1, x2, x3 + 1);
    f[iMM0] = (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3);
    f[iPM0] = (*this->splitB)(ePM0, x1, x2 + 1, x3);
    f[iM0M] = (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1);
    f[iP0M] = (*this->splitB)(eP0M, x1, x2, x3 + 1);
    f[i0MM] = (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1);
    f[i0PM] = (*this->splitB)(e0PM, x1, x2, x3 + 1);
    f[iMMM] = (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1);
    f[iPMM] = (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1);
    f[iMPM] = (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1);
    f[iPPM] = (*this->splitB)(ePPM, x1, x2, x3 + 1);

    f[d000] = (*this->split0)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    (*this->splitA)(eP00, x1, x2, x3) = f[dP00];
    (*this->splitA)(e0P0, x1, x2, x3) = f[d0P0];
    (*this->splitA)(e00P, x1, x2, x3) = f[d00P];
    (*this->splitA)(ePP0, x1, x2, x3) = f[dPP0];
    (*this->splitA)(eMP0, x1 + 1, x2, x3) = f[dMP0];
    (*this->splitA)(eP0P, x1, x2, x3) = f[dP0P];
    (*this->splitA)(eM0P, x1 + 1, x2, x3) = f[dM0P];
    (*this->splitA)(e0PP, x1, x2, x3) = f[d0PP];
    (*this->splitA)(e0MP, x1, x2 + 1, x3) = f[d0MP];
    (*this->splitA)(ePPP, x1, x2, x3) = f[dPPP];
    (*this->splitA)(eMPP, x1 + 1, x2, x3) = f[dMPP];
    (*this->splitA)(ePMP, x1, x2 + 1, x3) = f[dPMP];
    (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = f[dMMP];

    (*this->splitB)(eM00, x1 + 1, x2, x3) = f[dM00];
    (*this->splitB)(e0M0, x1, x2 + 1, x3) = f[d0M0];
    (*this->splitB)(e00M, x1, x2, x3 + 1) = f[d00M];
    (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = f[dMM0];
    (*this->splitB)(ePM0, x1, x2 + 1, x3) = f[dPM0];
    (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = f[dM0M];
    (*this->splitB)(eP0M, x1, x2, x3 + 1) = f[dP0M];
    (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = f[d0MM];
    (*this->splitB)(e0PM, x1, x2, x3 + 1) = f[d0PM];
    (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[dMMM];
    (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = f[dPMM];
    (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = f[dMPM];
    (*this->splitB)(ePPM, x1, x2, x3 + 1) = f[dPPM];

    (*this->split0)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                                unsigned long int direction)
{
    using namespace vf::lbm::dir;

    if ((direction & etP00) == etP00)
        (*this->splitB)(eM00, x1 + 1, x2, x3) = f[dP00];
    if ((direction & etM00) == etM00)
        (*this->splitA)(eP00, x1, x2, x3) = f[dM00];
    if ((direction & et0M0) == et0M0)
        (*this->splitA)(e0P0, x1, x2, x3) = f[d0M0];
    if ((direction & et0P0) == et0P0)
        (*this->splitB)(e0M0, x1, x2 + 1, x3) = f[d0P0];
    if ((direction & et00M) == et00M)
        (*this->splitA)(e00P, x1, x2, x3) = f[d00M];
    if ((direction & et00P) == et00P)
        (*this->splitB)(e00M, x1, x2, x3 + 1) = f[d00P];
    if ((direction & etMM0) == etMM0)
        (*this->splitA)(ePP0, x1, x2, x3) = f[dMM0];
    if ((direction & etPP0) == etPP0)
        (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = f[dPP0];
    if ((direction & etMP0) == etMP0)
        (*this->splitB)(ePM0, x1, x2 + 1, x3) = f[dMP0];
    if ((direction & etPM0) == etPM0)
        (*this->splitA)(eMP0, x1 + 1, x2, x3) = f[dPM0];
    if ((direction & etM0M) == etM0M)
        (*this->splitA)(eP0P, x1, x2, x3) = f[dM0M];
    if ((direction & etP0P) == etP0P)
        (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = f[dP0P];
    if ((direction & etM0P) == etM0P)
        (*this->splitB)(eP0M, x1, x2, x3 + 1) = f[dM0P];
    if ((direction & etP0M) == etP0M)
        (*this->splitA)(eM0P, x1 + 1, x2, x3) = f[dP0M];
    if ((direction & et0MM) == et0MM)
        (*this->splitA)(e0PP, x1, x2, x3) = f[d0MM];
    if ((direction & et0PP) == et0PP)
        (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = f[d0PP];
    if ((direction & et0MP) == et0MP)
        (*this->splitB)(e0PM, x1, x2, x3 + 1) = f[d0MP];
    if ((direction & et0PM) == et0PM)
        (*this->splitA)(e0MP, x1, x2 + 1, x3) = f[d0PM];
    if ((direction & etMMM) == etMMM)
        (*this->splitA)(ePPP, x1, x2, x3) = f[dMMM];
    if ((direction & etPPP) == etPPP)
        (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[dPPP];
    if ((direction & etPMM) == etPMM)
        (*this->splitA)(eMPP, x1 + 1, x2, x3) = f[dPMM];
    if ((direction & etMPP) == etMPP)
        (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = f[dMPP];
    if ((direction & etMPM) == etMPM)
        (*this->splitA)(ePMP, x1, x2 + 1, x3) = f[dMPM];
    if ((direction & etPMP) == etPMP)
        (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = f[dPMP];
    if ((direction & etPPM) == etPPM)
        (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = f[dPPM];
    if ((direction & etMMP) == etMMP)
        (*this->splitB)(ePPM, x1, x2, x3 + 1) = f[dMMP];
    if ((direction & et000) == et000)
        (*this->split0)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                                int direction)
{
    using namespace vf::lbm::dir;
 
    switch (direction) {
        case dP00:
            (*this->splitB)(eM00, x1 + 1, x2, x3) = f;
            break;
        case dM00:
            (*this->splitA)(eP00, x1, x2, x3) = f;
            break;
        case d0M0:
            (*this->splitA)(e0P0, x1, x2, x3) = f;
            break;
        case d0P0:
            (*this->splitB)(e0M0, x1, x2 + 1, x3) = f;
            break;
        case d00M:
            (*this->splitA)(e00P, x1, x2, x3) = f;
            break;
        case d00P:
            (*this->splitB)(e00M, x1, x2, x3 + 1) = f;
            break;
        case dMM0:
            (*this->splitA)(ePP0, x1, x2, x3) = f;
            break;
        case dPP0:
            (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = f;
            break;
        case dMP0:
            (*this->splitB)(ePM0, x1, x2 + 1, x3) = f;
            break;
        case dPM0:
            (*this->splitA)(eMP0, x1 + 1, x2, x3) = f;
            break;
        case dM0M:
            (*this->splitA)(eP0P, x1, x2, x3) = f;
            break;
        case dP0P:
            (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = f;
            break;
        case dM0P:
            (*this->splitB)(eP0M, x1, x2, x3 + 1) = f;
            break;
        case dP0M:
            (*this->splitA)(eM0P, x1 + 1, x2, x3) = f;
            break;
        case d0MM:
            (*this->splitA)(e0PP, x1, x2, x3) = f;
            break;
        case d0PP:
            (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = f;
            break;
        case d0MP:
            (*this->splitB)(e0PM, x1, x2, x3 + 1) = f;
            break;
        case d0PM:
            (*this->splitA)(e0MP, x1, x2 + 1, x3) = f;
            break;
        case dMMM:
            (*this->splitA)(ePPP, x1, x2, x3) = f;
            break;
        case dPPP:
            (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f;
            break;
        case dPMM:
            (*this->splitA)(eMPP, x1 + 1, x2, x3) = f;
            break;
        case dMPP:
            (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = f;
            break;
        case dMPM:
            (*this->splitA)(ePMP, x1, x2 + 1, x3) = f;
            break;
        case dPMP:
            (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = f;
            break;
        case dPPM:
            (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = f;
            break;
        case dMMP:
            (*this->splitB)(ePPM, x1, x2, x3 + 1) = f;
            break;
        case d000:
            (*this->split0)(x1, x2, x3) = f;
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
        (*this->splitA)(eP00, x1, x2, x3) = f[dP00];
    if ((direction & etM00) == etM00)
        (*this->splitB)(eM00, x1 + 1, x2, x3) = f[dM00];
    if ((direction & et0M0) == et0M0)
        (*this->splitB)(e0M0, x1, x2 + 1, x3) = f[d0M0];
    if ((direction & et0P0) == et0P0)
        (*this->splitA)(e0P0, x1, x2, x3) = f[d0P0];
    if ((direction & et00M) == et00M)
        (*this->splitB)(e00M, x1, x2, x3 + 1) = f[d00M];
    if ((direction & et00P) == et00P)
        (*this->splitA)(e00P, x1, x2, x3) = f[d00P];
    if ((direction & etMM0) == etMM0)
        (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = f[dMM0];
    if ((direction & etPP0) == etPP0)
        (*this->splitA)(ePP0, x1, x2, x3) = f[dPP0];
    if ((direction & etMP0) == etMP0)
        (*this->splitA)(eMP0, x1 + 1, x2, x3) = f[dMP0];
    if ((direction & etPM0) == etPM0)
        (*this->splitB)(ePM0, x1, x2 + 1, x3) = f[dPM0];
    if ((direction & etM0M) == etM0M)
        (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = f[dM0M];
    if ((direction & etP0P) == etP0P)
        (*this->splitA)(eP0P, x1, x2, x3) = f[dP0P];
    if ((direction & etM0P) == etM0P)
        (*this->splitA)(eM0P, x1 + 1, x2, x3) = f[dM0P];
    if ((direction & etP0M) == etP0M)
        (*this->splitB)(eP0M, x1, x2, x3 + 1) = f[dP0M];
    if ((direction & et0MM) == et0MM)
        (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = f[d0MM];
    if ((direction & et0PP) == et0PP)
        (*this->splitA)(e0PP, x1, x2, x3) = f[d0PP];
    if ((direction & et0MP) == et0MP)
        (*this->splitA)(e0MP, x1, x2 + 1, x3) = f[d0MP];
    if ((direction & et0PM) == et0PM)
        (*this->splitB)(e0PM, x1, x2, x3 + 1) = f[d0PM];
    if ((direction & etMMM) == etMMM)
        (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f[dMMM];
    if ((direction & etPPP) == etPPP)
        (*this->splitA)(ePPP, x1, x2, x3) = f[dPPP];
    if ((direction & etPMM) == etPMM)
        (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = f[dPMM];
    if ((direction & etMPP) == etMPP)
        (*this->splitA)(eMPP, x1 + 1, x2, x3) = f[dMPP];
    if ((direction & etMPM) == etMPM)
        (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = f[dMPM];
    if ((direction & etPMP) == etPMP)
        (*this->splitA)(ePMP, x1, x2 + 1, x3) = f[dPMP];
    if ((direction & etPPM) == etPPM)
        (*this->splitB)(ePPM, x1, x2, x3 + 1) = f[dPPM];
    if ((direction & etMMP) == etMMP)
        (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = f[dMMP];
    if ((direction & et000) == et000)
        (*this->split0)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                                   unsigned long int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dP00:
            (*this->splitA)(eP00, x1, x2, x3) = f;
            break;
        case dM00:
            (*this->splitB)(eM00, x1 + 1, x2, x3) = f;
            break;
        case d0M0:
            (*this->splitB)(e0M0, x1, x2 + 1, x3) = f;
            break;
        case d0P0:
            (*this->splitA)(e0P0, x1, x2, x3) = f;
            break;
        case d00M:
            (*this->splitB)(e00M, x1, x2, x3 + 1) = f;
            break;
        case d00P:
            (*this->splitA)(e00P, x1, x2, x3) = f;
            break;
        case dMM0:
            (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = f;
            break;
        case dPP0:
            (*this->splitA)(ePP0, x1, x2, x3) = f;
            break;
        case dMP0:
            (*this->splitA)(eMP0, x1 + 1, x2, x3) = f;
            break;
        case dPM0:
            (*this->splitB)(ePM0, x1, x2 + 1, x3) = f;
            break;
        case dM0M:
            (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = f;
            break;
        case dP0P:
            (*this->splitA)(eP0P, x1, x2, x3) = f;
            break;
        case dM0P:
            (*this->splitA)(eM0P, x1 + 1, x2, x3) = f;
            break;
        case dP0M:
            (*this->splitB)(eP0M, x1, x2, x3 + 1) = f;
            break;
        case d0MM:
            (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = f;
            break;
        case d0PP:
            (*this->splitA)(e0PP, x1, x2, x3) = f;
            break;
        case d0MP:
            (*this->splitA)(e0MP, x1, x2 + 1, x3) = f;
            break;
        case d0PM:
            (*this->splitB)(e0PM, x1, x2, x3 + 1) = f;
            break;
        case dMMM:
            (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = f;
            break;
        case dPPP:
            (*this->splitA)(ePPP, x1, x2, x3) = f;
            break;
        case dPMM:
            (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = f;
            break;
        case dMPP:
            (*this->splitA)(eMPP, x1 + 1, x2, x3) = f;
            break;
        case dMPM:
            (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = f;
            break;
        case dPMP:
            (*this->splitA)(ePMP, x1, x2 + 1, x3) = f;
            break;
        case dPPM:
            (*this->splitB)(ePPM, x1, x2, x3 + 1) = f;
            break;
        case dMMP:
            (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = f;
            break;
        case d000:
            (*this->split0)(x1, x2, x3) = f;
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
            return (*this->splitB)(eM00, x1 + 1, x2, x3);
        case dP00:
            return (*this->splitA)(eP00, x1, x2, x3);
        case d0P0:
            return (*this->splitA)(e0P0, x1, x2, x3);
        case d0M0:
            return (*this->splitB)(e0M0, x1, x2 + 1, x3);
        case d00P:
            return (*this->splitA)(e00P, x1, x2, x3);
        case d00M:
            return (*this->splitB)(e00M, x1, x2, x3 + 1);
        case dPP0:
            return (*this->splitA)(ePP0, x1, x2, x3);
        case dMM0:
            return (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3);
        case dPM0:
            return (*this->splitB)(ePM0, x1, x2 + 1, x3);
        case dMP0:
            return (*this->splitA)(eMP0, x1 + 1, x2, x3);
        case dP0P:
            return (*this->splitA)(eP0P, x1, x2, x3);
        case dM0M:
            return (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1);
        case dP0M:
            return (*this->splitB)(eP0M, x1, x2, x3 + 1);
        case dM0P:
            return (*this->splitA)(eM0P, x1 + 1, x2, x3);
        case d0PP:
            return (*this->splitA)(e0PP, x1, x2, x3);
        case d0MM:
            return (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1);
        case d0PM:
            return (*this->splitB)(e0PM, x1, x2, x3 + 1);
        case d0MP:
            return (*this->splitA)(e0MP, x1, x2 + 1, x3);
        case dPPP:
            return (*this->splitA)(ePPP, x1, x2, x3);
        case dMMM:
            return (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1);
        case dMPP:
            return (*this->splitA)(eMPP, x1 + 1, x2, x3);
        case dPMM:
            return (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1);
        case dPMP:
            return (*this->splitA)(ePMP, x1, x2 + 1, x3);
        case dMPM:
            return (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1);
        case dMMP:
            return (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3);
        case dPPM:
            return (*this->splitB)(ePPM, x1, x2, x3 + 1);
        case d000:
            return (*this->split0)(x1, x2, x3);
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
            return (*this->splitB)(eM00, x1 + 1, x2, x3);
        case dM00:
            return (*this->splitA)(eP00, x1, x2, x3);
        case d0M0:
            return (*this->splitA)(e0P0, x1, x2, x3);
        case d0P0:
            return (*this->splitB)(e0M0, x1, x2 + 1, x3);
        case d00M:
            return (*this->splitA)(e00P, x1, x2, x3);
        case d00P:
            return (*this->splitB)(e00M, x1, x2, x3 + 1);
        case dMM0:
            return (*this->splitA)(ePP0, x1, x2, x3);
        case dPP0:
            return (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3);
        case dMP0:
            return (*this->splitB)(ePM0, x1, x2 + 1, x3);
        case dPM0:
            return (*this->splitA)(eMP0, x1 + 1, x2, x3);
        case dM0M:
            return (*this->splitA)(eP0P, x1, x2, x3);
        case dP0P:
            return (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1);
        case dM0P:
            return (*this->splitB)(eP0M, x1, x2, x3 + 1);
        case dP0M:
            return (*this->splitA)(eM0P, x1 + 1, x2, x3);
        case d0MM:
            return (*this->splitA)(e0PP, x1, x2, x3);
        case d0PP:
            return (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1);
        case d0MP:
            return (*this->splitB)(e0PM, x1, x2, x3 + 1);
        case d0PM:
            return (*this->splitA)(e0MP, x1, x2 + 1, x3);
        case dMMM:
            return (*this->splitA)(ePPP, x1, x2, x3);
        case dPPP:
            return (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1);
        case dPMM:
            return (*this->splitA)(eMPP, x1 + 1, x2, x3);
        case dMPP:
            return (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1);
        case dMPM:
            return (*this->splitA)(ePMP, x1, x2 + 1, x3);
        case dPMP:
            return (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1);
        case dPPM:
            return (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3);
        case dMMP:
            return (*this->splitB)(ePPM, x1, x2, x3 + 1);
        case d000:
            return (*this->split0)(x1, x2, x3);
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
CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr EsoSplit::getSplitA()
{
    return this->splitA;
}
//////////////////////////////////////////////////////////////////////////
CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr EsoSplit::getSplitB()
{
    return this->splitB;
}
//////////////////////////////////////////////////////////////////////////
CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr EsoSplit::getSplit0()
{
    return this->split0;
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNX1(size_t newNX1) { NX1 = newNX1; }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNX2(size_t newNX2) { NX2 = newNX2; }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setNX3(size_t newNX3) { NX3 = newNX3; }
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setSplitA(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array)
{
    splitA = array;
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setSplitB(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array)
{
    splitB = array;
}
//////////////////////////////////////////////////////////////////////////
void EsoSplit::setSplit0(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr array)
{
    split0 = array;
}

//////////////////////////////////////////////////////////////////////////

//! \}
