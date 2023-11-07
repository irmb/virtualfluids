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
//! \file D3Q27EsoTwist3DSplittedVector.cpp
//! \ingroup Data
//! \author Konstantin Kutscher
//=======================================================================================

#include "D3Q27EsoTwist3DSplittedVector.h"
#include "EsoTwistD3Q27System.h"

D3Q27EsoTwist3DSplittedVector::D3Q27EsoTwist3DSplittedVector() = default;
//////////////////////////////////////////////////////////////////////////
D3Q27EsoTwist3DSplittedVector::D3Q27EsoTwist3DSplittedVector(size_t nx1, size_t nx2, size_t nx3, real value)
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
D3Q27EsoTwist3DSplittedVector::~D3Q27EsoTwist3DSplittedVector() = default;
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::swap() { std::swap(this->localDistributions, this->nonLocalDistributions); }
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    const size_t x1p = x1 + 1;
    const size_t x2p = x2 + 1;
    const size_t x3p = x3 + 1;

    f[vf::lbm::dir::dP00] = (*this->localDistributions)(D3Q27System::ET_P00, x1, x2, x3);
    f[vf::lbm::dir::DIR_0P0] = (*this->localDistributions)(D3Q27System::ET_0P0, x1, x2, x3);
    f[vf::lbm::dir::DIR_00P] = (*this->localDistributions)(D3Q27System::ET_00P, x1, x2, x3);
    f[vf::lbm::dir::DIR_PP0] = (*this->localDistributions)(D3Q27System::ET_PP0, x1, x2, x3);
    f[vf::lbm::dir::DIR_MP0] = (*this->localDistributions)(D3Q27System::ET_MP0, x1p, x2, x3);
    f[vf::lbm::dir::DIR_P0P] = (*this->localDistributions)(D3Q27System::ET_P0P, x1, x2, x3);
    f[vf::lbm::dir::DIR_M0P] = (*this->localDistributions)(D3Q27System::ET_M0P, x1p, x2, x3);
    f[vf::lbm::dir::DIR_0PP] = (*this->localDistributions)(D3Q27System::ET_0PP, x1, x2, x3);
    f[vf::lbm::dir::DIR_0MP] = (*this->localDistributions)(D3Q27System::ET_0MP, x1, x2p, x3);
    f[vf::lbm::dir::DIR_PPP] = (*this->localDistributions)(D3Q27System::ET_PPP, x1, x2, x3);
    f[vf::lbm::dir::DIR_MPP] = (*this->localDistributions)(D3Q27System::ET_MPP, x1p, x2, x3);
    f[vf::lbm::dir::DIR_PMP] = (*this->localDistributions)(D3Q27System::ET_PMP, x1, x2p, x3);
    f[vf::lbm::dir::DIR_MMP] = (*this->localDistributions)(D3Q27System::ET_MMP, x1p, x2p, x3);

    f[vf::lbm::dir::dM00] = (*this->nonLocalDistributions)(D3Q27System::ET_M00, x1p, x2, x3);
    f[vf::lbm::dir::DIR_0M0] = (*this->nonLocalDistributions)(D3Q27System::ET_0M0, x1, x2p, x3);
    f[vf::lbm::dir::DIR_00M] = (*this->nonLocalDistributions)(D3Q27System::ET_00M, x1, x2, x3p);
    f[vf::lbm::dir::DIR_MM0] = (*this->nonLocalDistributions)(D3Q27System::ET_MM0, x1p, x2p, x3);
    f[vf::lbm::dir::DIR_PM0] = (*this->nonLocalDistributions)(D3Q27System::ET_PM0, x1, x2p, x3);
    f[vf::lbm::dir::DIR_M0M] = (*this->nonLocalDistributions)(D3Q27System::ET_M0M, x1p, x2, x3p);
    f[vf::lbm::dir::DIR_P0M] = (*this->nonLocalDistributions)(D3Q27System::ET_P0M, x1, x2, x3p);
    f[vf::lbm::dir::DIR_0MM] = (*this->nonLocalDistributions)(D3Q27System::ET_0MM, x1, x2p, x3p);
    f[vf::lbm::dir::DIR_0PM] = (*this->nonLocalDistributions)(D3Q27System::ET_0PM, x1, x2, x3p);
    f[vf::lbm::dir::DIR_MMM] = (*this->nonLocalDistributions)(D3Q27System::ET_MMM, x1p, x2p, x3p);
    f[vf::lbm::dir::DIR_PMM] = (*this->nonLocalDistributions)(D3Q27System::ET_PMM, x1, x2p, x3p);
    f[vf::lbm::dir::DIR_MPM] = (*this->nonLocalDistributions)(D3Q27System::ET_MPM, x1p, x2, x3p);
    f[vf::lbm::dir::DIR_PPM] = (*this->nonLocalDistributions)(D3Q27System::ET_PPM, x1, x2, x3p);

    f[vf::lbm::dir::d000] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    (*this->localDistributions)(D3Q27System::ET_P00, x1, x2, x3) = f[INV_P00];
    (*this->localDistributions)(D3Q27System::ET_0P0, x1, x2, x3) = f[INV_0P0];
    (*this->localDistributions)(D3Q27System::ET_00P, x1, x2, x3) = f[INV_00P];
    (*this->localDistributions)(D3Q27System::ET_PP0, x1, x2, x3) = f[INV_PP0];
    (*this->localDistributions)(D3Q27System::ET_MP0, x1 + 1, x2, x3) = f[INV_MP0];
    (*this->localDistributions)(D3Q27System::ET_P0P, x1, x2, x3) = f[INV_P0P];
    (*this->localDistributions)(D3Q27System::ET_M0P, x1 + 1, x2, x3) = f[INV_M0P];
    (*this->localDistributions)(D3Q27System::ET_0PP, x1, x2, x3) = f[INV_0PP];
    (*this->localDistributions)(D3Q27System::ET_0MP, x1, x2 + 1, x3) = f[INV_0MP];
    (*this->localDistributions)(D3Q27System::ET_PPP, x1, x2, x3) = f[INV_PPP];
    (*this->localDistributions)(D3Q27System::ET_MPP, x1 + 1, x2, x3) = f[INV_MPP];
    (*this->localDistributions)(D3Q27System::ET_PMP, x1, x2 + 1, x3) = f[INV_PMP];
    (*this->localDistributions)(D3Q27System::ET_MMP, x1 + 1, x2 + 1, x3) = f[INV_MMP];

    (*this->nonLocalDistributions)(D3Q27System::ET_M00, x1 + 1, x2, x3) = f[INV_M00];
    (*this->nonLocalDistributions)(D3Q27System::ET_0M0, x1, x2 + 1, x3) = f[INV_0M0];
    (*this->nonLocalDistributions)(D3Q27System::ET_00M, x1, x2, x3 + 1) = f[INV_00M];
    (*this->nonLocalDistributions)(D3Q27System::ET_MM0, x1 + 1, x2 + 1, x3) = f[INV_MM0];
    (*this->nonLocalDistributions)(D3Q27System::ET_PM0, x1, x2 + 1, x3) = f[INV_PM0];
    (*this->nonLocalDistributions)(D3Q27System::ET_M0M, x1 + 1, x2, x3 + 1) = f[INV_M0M];
    (*this->nonLocalDistributions)(D3Q27System::ET_P0M, x1, x2, x3 + 1) = f[INV_P0M];
    (*this->nonLocalDistributions)(D3Q27System::ET_0MM, x1, x2 + 1, x3 + 1) = f[INV_0MM];
    (*this->nonLocalDistributions)(D3Q27System::ET_0PM, x1, x2, x3 + 1) = f[INV_0PM];
    (*this->nonLocalDistributions)(D3Q27System::ET_MMM, x1 + 1, x2 + 1, x3 + 1) = f[INV_MMM];
    (*this->nonLocalDistributions)(D3Q27System::ET_PMM, x1, x2 + 1, x3 + 1) = f[INV_PMM];
    (*this->nonLocalDistributions)(D3Q27System::ET_MPM, x1 + 1, x2, x3 + 1) = f[INV_MPM];
    (*this->nonLocalDistributions)(D3Q27System::ET_PPM, x1, x2, x3 + 1) = f[INV_PPM];

    (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    f[INV_P00] = (*this->localDistributions)(D3Q27System::ET_P00, x1, x2, x3);
    f[INV_0P0] = (*this->localDistributions)(D3Q27System::ET_0P0, x1, x2, x3);
    f[INV_00P] = (*this->localDistributions)(D3Q27System::ET_00P, x1, x2, x3);
    f[INV_PP0] = (*this->localDistributions)(D3Q27System::ET_PP0, x1, x2, x3);
    f[INV_MP0] = (*this->localDistributions)(D3Q27System::ET_MP0, x1 + 1, x2, x3);
    f[INV_P0P] = (*this->localDistributions)(D3Q27System::ET_P0P, x1, x2, x3);
    f[INV_M0P] = (*this->localDistributions)(D3Q27System::ET_M0P, x1 + 1, x2, x3);
    f[INV_0PP] = (*this->localDistributions)(D3Q27System::ET_0PP, x1, x2, x3);
    f[INV_0MP] = (*this->localDistributions)(D3Q27System::ET_0MP, x1, x2 + 1, x3);
    f[INV_PPP] = (*this->localDistributions)(D3Q27System::ET_PPP, x1, x2, x3);
    f[INV_MPP] = (*this->localDistributions)(D3Q27System::ET_MPP, x1 + 1, x2, x3);
    f[INV_PMP] = (*this->localDistributions)(D3Q27System::ET_PMP, x1, x2 + 1, x3);
    f[INV_MMP] = (*this->localDistributions)(D3Q27System::ET_MMP, x1 + 1, x2 + 1, x3);

    f[INV_M00] = (*this->nonLocalDistributions)(D3Q27System::ET_M00, x1 + 1, x2, x3);
    f[INV_0M0] = (*this->nonLocalDistributions)(D3Q27System::ET_0M0, x1, x2 + 1, x3);
    f[INV_00M] = (*this->nonLocalDistributions)(D3Q27System::ET_00M, x1, x2, x3 + 1);
    f[INV_MM0] = (*this->nonLocalDistributions)(D3Q27System::ET_MM0, x1 + 1, x2 + 1, x3);
    f[INV_PM0] = (*this->nonLocalDistributions)(D3Q27System::ET_PM0, x1, x2 + 1, x3);
    f[INV_M0M] = (*this->nonLocalDistributions)(D3Q27System::ET_M0M, x1 + 1, x2, x3 + 1);
    f[INV_P0M] = (*this->nonLocalDistributions)(D3Q27System::ET_P0M, x1, x2, x3 + 1);
    f[INV_0MM] = (*this->nonLocalDistributions)(D3Q27System::ET_0MM, x1, x2 + 1, x3 + 1);
    f[INV_0PM] = (*this->nonLocalDistributions)(D3Q27System::ET_0PM, x1, x2, x3 + 1);
    f[INV_MMM] = (*this->nonLocalDistributions)(D3Q27System::ET_MMM, x1 + 1, x2 + 1, x3 + 1);
    f[INV_PMM] = (*this->nonLocalDistributions)(D3Q27System::ET_PMM, x1, x2 + 1, x3 + 1);
    f[INV_MPM] = (*this->nonLocalDistributions)(D3Q27System::ET_MPM, x1 + 1, x2, x3 + 1);
    f[INV_PPM] = (*this->nonLocalDistributions)(D3Q27System::ET_PPM, x1, x2, x3 + 1);

    f[d000] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    (*this->localDistributions)(D3Q27System::ET_P00, x1, x2, x3) = f[dP00];
    (*this->localDistributions)(D3Q27System::ET_0P0, x1, x2, x3) = f[DIR_0P0];
    (*this->localDistributions)(D3Q27System::ET_00P, x1, x2, x3) = f[DIR_00P];
    (*this->localDistributions)(D3Q27System::ET_PP0, x1, x2, x3) = f[DIR_PP0];
    (*this->localDistributions)(D3Q27System::ET_MP0, x1 + 1, x2, x3) = f[DIR_MP0];
    (*this->localDistributions)(D3Q27System::ET_P0P, x1, x2, x3) = f[DIR_P0P];
    (*this->localDistributions)(D3Q27System::ET_M0P, x1 + 1, x2, x3) = f[DIR_M0P];
    (*this->localDistributions)(D3Q27System::ET_0PP, x1, x2, x3) = f[DIR_0PP];
    (*this->localDistributions)(D3Q27System::ET_0MP, x1, x2 + 1, x3) = f[DIR_0MP];
    (*this->localDistributions)(D3Q27System::ET_PPP, x1, x2, x3) = f[DIR_PPP];
    (*this->localDistributions)(D3Q27System::ET_MPP, x1 + 1, x2, x3) = f[DIR_MPP];
    (*this->localDistributions)(D3Q27System::ET_PMP, x1, x2 + 1, x3) = f[DIR_PMP];
    (*this->localDistributions)(D3Q27System::ET_MMP, x1 + 1, x2 + 1, x3) = f[DIR_MMP];

    (*this->nonLocalDistributions)(D3Q27System::ET_M00, x1 + 1, x2, x3) = f[dM00];
    (*this->nonLocalDistributions)(D3Q27System::ET_0M0, x1, x2 + 1, x3) = f[DIR_0M0];
    (*this->nonLocalDistributions)(D3Q27System::ET_00M, x1, x2, x3 + 1) = f[DIR_00M];
    (*this->nonLocalDistributions)(D3Q27System::ET_MM0, x1 + 1, x2 + 1, x3) = f[DIR_MM0];
    (*this->nonLocalDistributions)(D3Q27System::ET_PM0, x1, x2 + 1, x3) = f[DIR_PM0];
    (*this->nonLocalDistributions)(D3Q27System::ET_M0M, x1 + 1, x2, x3 + 1) = f[DIR_M0M];
    (*this->nonLocalDistributions)(D3Q27System::ET_P0M, x1, x2, x3 + 1) = f[DIR_P0M];
    (*this->nonLocalDistributions)(D3Q27System::ET_0MM, x1, x2 + 1, x3 + 1) = f[DIR_0MM];
    (*this->nonLocalDistributions)(D3Q27System::ET_0PM, x1, x2, x3 + 1) = f[DIR_0PM];
    (*this->nonLocalDistributions)(D3Q27System::ET_MMM, x1 + 1, x2 + 1, x3 + 1) = f[DIR_MMM];
    (*this->nonLocalDistributions)(D3Q27System::ET_PMM, x1, x2 + 1, x3 + 1) = f[DIR_PMM];
    (*this->nonLocalDistributions)(D3Q27System::ET_MPM, x1 + 1, x2, x3 + 1) = f[DIR_MPM];
    (*this->nonLocalDistributions)(D3Q27System::ET_PPM, x1, x2, x3 + 1) = f[DIR_PPM];

    (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                                unsigned long int direction)
{
    using namespace vf::lbm::dir;

    if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
        (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = f[dP00];
    if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
        (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = f[dM00];
    if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
        (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = f[DIR_0M0];
    if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
        (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = f[DIR_0P0];
    if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
        (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = f[DIR_00M];
    if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
        (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = f[DIR_00P];
    if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
        (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = f[DIR_MM0];
    if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
        (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = f[DIR_PP0];
    if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
        (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = f[DIR_MP0];
    if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
        (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = f[DIR_PM0];
    if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
        (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = f[DIR_M0M];
    if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
        (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = f[DIR_P0P];
    if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
        (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = f[DIR_M0P];
    if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
        (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = f[DIR_P0M];
    if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
        (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = f[DIR_0MM];
    if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
        (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = f[DIR_0PP];
    if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
        (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = f[DIR_0MP];
    if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
        (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = f[DIR_0PM];
    if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
        (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = f[DIR_MMM];
    if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
        (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = f[DIR_PPP];
    if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
        (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = f[DIR_PMM];
    if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
        (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = f[DIR_MPP];
    if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
        (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = f[DIR_MPM];
    if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
        (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = f[DIR_PMP];
    if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
        (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = f[DIR_PPM];
    if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
        (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = f[DIR_MMP];
    if ((direction & EsoTwistD3Q27System::REST) == EsoTwistD3Q27System::REST)
        (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                                int direction)
{
    using namespace vf::lbm::dir;
 
    switch (direction) {
        case dP00:
            (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = f;
            break;
        case dM00:
            (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = f;
            break;
        case DIR_0M0:
            (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = f;
            break;
        case DIR_0P0:
            (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = f;
            break;
        case DIR_00M:
            (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = f;
            break;
        case DIR_00P:
            (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = f;
            break;
        case DIR_MM0:
            (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = f;
            break;
        case DIR_PP0:
            (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = f;
            break;
        case DIR_MP0:
            (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = f;
            break;
        case DIR_PM0:
            (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = f;
            break;
        case DIR_M0M:
            (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = f;
            break;
        case DIR_P0P:
            (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = f;
            break;
        case DIR_M0P:
            (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = f;
            break;
        case DIR_P0M:
            (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = f;
            break;
        case DIR_0MM:
            (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = f;
            break;
        case DIR_0PP:
            (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = f;
            break;
        case DIR_0MP:
            (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = f;
            break;
        case DIR_0PM:
            (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = f;
            break;
        case DIR_MMM:
            (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = f;
            break;
        case DIR_PPP:
            (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = f;
            break;
        case DIR_PMM:
            (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = f;
            break;
        case DIR_MPP:
            (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = f;
            break;
        case DIR_MPM:
            (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = f;
            break;
        case DIR_PMP:
            (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = f;
            break;
        case DIR_PPM:
            (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = f;
            break;
        case DIR_MMP:
            (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = f;
            break;
        case d000:
            (*this->zeroDistributions)(x1, x2, x3) = f;
            break;
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2,
                                                                   size_t x3, unsigned long int direction)
{
    using namespace vf::lbm::dir;

    if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
        (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = f[dP00];
    if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
        (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = f[dM00];
    if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
        (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = f[DIR_0M0];
    if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
        (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = f[DIR_0P0];
    if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
        (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = f[DIR_00M];
    if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
        (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = f[DIR_00P];
    if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
        (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = f[DIR_MM0];
    if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
        (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = f[DIR_PP0];
    if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
        (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = f[DIR_MP0];
    if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
        (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = f[DIR_PM0];
    if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
        (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = f[DIR_M0M];
    if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
        (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = f[DIR_P0P];
    if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
        (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = f[DIR_M0P];
    if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
        (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = f[DIR_P0M];
    if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
        (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = f[DIR_0MM];
    if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
        (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = f[DIR_0PP];
    if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
        (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = f[DIR_0MP];
    if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
        (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = f[DIR_0PM];
    if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
        (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = f[DIR_MMM];
    if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
        (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = f[DIR_PPP];
    if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
        (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = f[DIR_PMM];
    if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
        (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = f[DIR_MPP];
    if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
        (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = f[DIR_MPM];
    if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
        (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = f[DIR_PMP];
    if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
        (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = f[DIR_PPM];
    if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
        (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = f[DIR_MMP];
    if ((direction & EsoTwistD3Q27System::REST) == EsoTwistD3Q27System::REST)
        (*this->zeroDistributions)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                                   unsigned long int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dP00:
            (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = f;
            break;
        case dM00:
            (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = f;
            break;
        case DIR_0M0:
            (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = f;
            break;
        case DIR_0P0:
            (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = f;
            break;
        case DIR_00M:
            (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = f;
            break;
        case DIR_00P:
            (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = f;
            break;
        case DIR_MM0:
            (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = f;
            break;
        case DIR_PP0:
            (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = f;
            break;
        case DIR_MP0:
            (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = f;
            break;
        case DIR_PM0:
            (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = f;
            break;
        case DIR_M0M:
            (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = f;
            break;
        case DIR_P0P:
            (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = f;
            break;
        case DIR_M0P:
            (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = f;
            break;
        case DIR_P0M:
            (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = f;
            break;
        case DIR_0MM:
            (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = f;
            break;
        case DIR_0PP:
            (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = f;
            break;
        case DIR_0MP:
            (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = f;
            break;
        case DIR_0PM:
            (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = f;
            break;
        case DIR_MMM:
            (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = f;
            break;
        case DIR_PPP:
            (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = f;
            break;
        case DIR_PMM:
            (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = f;
            break;
        case DIR_MPP:
            (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = f;
            break;
        case DIR_MPM:
            (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = f;
            break;
        case DIR_PMP:
            (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = f;
            break;
        case DIR_PPM:
            (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = f;
            break;
        case DIR_MMP:
            (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = f;
            break;
        case d000:
            (*this->zeroDistributions)(x1, x2, x3) = f;
            break;
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
real D3Q27EsoTwist3DSplittedVector::getPreCollisionDistributionForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dM00:
            return (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
        case dP00:
            return (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
        case DIR_0P0:
            return (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
        case DIR_0M0:
            return (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
        case DIR_00P:
            return (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
        case DIR_00M:
            return (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
        case DIR_PP0:
            return (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
        case DIR_MM0:
            return (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
        case DIR_PM0:
            return (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
        case DIR_MP0:
            return (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
        case DIR_P0P:
            return (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
        case DIR_M0M:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
        case DIR_P0M:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
        case DIR_M0P:
            return (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
        case DIR_0PP:
            return (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
        case DIR_0MM:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
        case DIR_0PM:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
        case DIR_0MP:
            return (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
        case DIR_PPP:
            return (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
        case DIR_MMM:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
        case DIR_MPP:
            return (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
        case DIR_PMM:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
        case DIR_PMP:
            return (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
        case DIR_MPM:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
        case DIR_MMP:
            return (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);
        case DIR_PPM:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);
        case d000:
            return (*this->zeroDistributions)(x1, x2, x3);
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
real D3Q27EsoTwist3DSplittedVector::getDistributionInvForDirection(size_t x1, size_t x2, size_t x3, int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dP00:
            return (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
        case dM00:
            return (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
        case DIR_0M0:
            return (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
        case DIR_0P0:
            return (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
        case DIR_00M:
            return (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
        case DIR_00P:
            return (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
        case DIR_MM0:
            return (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
        case DIR_PP0:
            return (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
        case DIR_MP0:
            return (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
        case DIR_PM0:
            return (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
        case DIR_M0M:
            return (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
        case DIR_P0P:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
        case DIR_M0P:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
        case DIR_P0M:
            return (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
        case DIR_0MM:
            return (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
        case DIR_0PP:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
        case DIR_0MP:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
        case DIR_0PM:
            return (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
        case DIR_MMM:
            return (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
        case DIR_PPP:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
        case DIR_PMM:
            return (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
        case DIR_MPP:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
        case DIR_MPM:
            return (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
        case DIR_PMP:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
        case DIR_PPM:
            return (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);
        case DIR_MMP:
            return (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);
        case d000:
            return (*this->zeroDistributions)(x1, x2, x3);
        default:
            UB_THROW(UbException(UB_EXARGS, "Direction didn't find"));
    }
}
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSplittedVector::getNX1() const { return NX1; }
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSplittedVector::getNX2() const { return NX2; }
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSplittedVector::getNX3() const { return NX3; }
//////////////////////////////////////////////////////////////////////////
CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr D3Q27EsoTwist3DSplittedVector::getLocalDistributions()
{
    return this->localDistributions;
}
//////////////////////////////////////////////////////////////////////////
CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr D3Q27EsoTwist3DSplittedVector::getNonLocalDistributions()
{
    return this->nonLocalDistributions;
}
//////////////////////////////////////////////////////////////////////////
CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr D3Q27EsoTwist3DSplittedVector::getZeroDistributions()
{
    return this->zeroDistributions;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setNX1(size_t newNX1) { NX1 = newNX1; }
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setNX2(size_t newNX2) { NX2 = newNX2; }
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setNX3(size_t newNX3) { NX3 = newNX3; }
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array)
{
    localDistributions = array;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setNonLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr array)
{
    nonLocalDistributions = array;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSplittedVector::setZeroDistributions(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr array)
{
    zeroDistributions = array;
}

//////////////////////////////////////////////////////////////////////////
