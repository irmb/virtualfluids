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
//! \file OneDistributionFullDirectConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef OneDistributionFullDirectConnector_H
#define OneDistributionFullDirectConnector_H

#include "Block3D.h"
#include "D3Q27System.h"
#include "FullDirectConnector.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

//! \brief   Exchange data between blocks.
//! \details Connector send and receive full distributions between two blocks in shared memory.
class OneDistributionFullDirectConnector : public FullDirectConnector
{
public:
    OneDistributionFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);
    void init() override;

protected:
    inline void updatePointers();
    inline void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);

private:
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFrom;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFrom;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsFrom;

    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTo;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTo;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsTo;

    SPtr<EsoTwist3D> fFrom;
    SPtr<EsoTwist3D> fTo;
};
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullDirectConnector::updatePointers()
{
    localDistributionsFrom = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getLocalDistributions();
    nonLocalDistributionsFrom = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getNonLocalDistributions();
    zeroDistributionsFrom = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getZeroDistributions();

    localDistributionsTo = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getLocalDistributions();
    nonLocalDistributionsTo = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getNonLocalDistributions();
    zeroDistributionsTo     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{
    (*this->localDistributionsTo)(D3Q27System::ET_E, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_E, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_N, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_N, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_T, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_T, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_NE, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_NE, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TE, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TE, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TN, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TN, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TNE, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TNE, x1From, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
    (*this->localDistributionsTo)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) =
        (*this->localDistributionsFrom)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

    (*this->nonLocalDistributionsTo)(D3Q27System::ET_W, x1To + 1, x2To, x3To) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_S, x1To, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_B, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

    (*this->zeroDistributionsTo)(x1To, x2To, x3To) = (*this->zeroDistributionsFrom)(x1From, x2From, x3From);
}
#endif
