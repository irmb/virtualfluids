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
//! \file ThreeDistributionsDoubleGhostLayerFullDirectConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef ThreeDistributionsDoubleGhostLayerFullDirectConnector_H
#define ThreeDistributionsDoubleGhostLayerFullDirectConnector_H

#include "FullDirectConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"
#include "DataSet3D.h"

//! \brief   Exchange data between blocks. 
//! \details Connector send and receive full distributions between two blocks in shared memory.

class ThreeDistributionsDoubleGhostLayerFullDirectConnector : public FullDirectConnector
{
public:
	ThreeDistributionsDoubleGhostLayerFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);
    void init() override;
    void sendVectors() override;

protected:
    inline void updatePointers() override;
    void exchangeData() override;
    inline void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To) override;

private:
	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromf;
	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromf;
	CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFromf;

	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTof;
	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTof;
	CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsTof;

	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromh;
	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromh;
	CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFromh;

	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsToh;
	CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsToh;
	CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsToh;

	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromh2;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromh2;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsFromh2;

    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsToh2;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsToh2;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsToh2;

	SPtr<EsoTwist3D> fFrom, hFrom, hFrom2;
    SPtr<EsoTwist3D> fTo, hTo, hTo2;

    SPtr<PressureFieldArray3D> pressureFrom, pressureTo;
};
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsDoubleGhostLayerFullDirectConnector::updatePointers()
{
    localDistributionsFromf = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getLocalDistributions();
    nonLocalDistributionsFromf = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getNonLocalDistributions();
    zeroDistributionsFromf = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getZeroDistributions();

    localDistributionsTof    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getLocalDistributions();
    nonLocalDistributionsTof = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getNonLocalDistributions();
    zeroDistributionsTof     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getZeroDistributions();

    localDistributionsFromh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom)->getLocalDistributions();
    nonLocalDistributionsFromh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom)->getNonLocalDistributions();
    zeroDistributionsFromh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom)->getZeroDistributions();

    localDistributionsToh    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo)->getLocalDistributions();
    nonLocalDistributionsToh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo)->getNonLocalDistributions();
    zeroDistributionsToh     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo)->getZeroDistributions();

    localDistributionsFromh2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom2)->getLocalDistributions();
    nonLocalDistributionsFromh2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom2)->getNonLocalDistributions();
    zeroDistributionsFromh2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom2)->getZeroDistributions();

    localDistributionsToh2    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo2)->getLocalDistributions();
    nonLocalDistributionsToh2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo2)->getNonLocalDistributions();
    zeroDistributionsToh2     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo2)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsDoubleGhostLayerFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{
	(*this->localDistributionsTof)(D3Q27System::ET_E, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_E, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_N, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_N, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_T, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_T, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_NE, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_NE, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TE, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TE, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TN, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TN, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TNE, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TNE, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

	(*this->nonLocalDistributionsTof)(D3Q27System::ET_W, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_S, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_B, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

	(*this->zeroDistributionsTof)(x1To, x2To, x3To) = (*this->zeroDistributionsFromf)(x1From, x2From, x3From);


	(*this->localDistributionsToh)(D3Q27System::ET_E, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_E, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_N, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_N, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_T, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_T, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_NE, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_NE, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TE, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TE, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TN, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TN, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TNE, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TNE, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

	(*this->nonLocalDistributionsToh)(D3Q27System::ET_W, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_S, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_B, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

	(*this->zeroDistributionsToh)(x1To, x2To, x3To) = (*this->zeroDistributionsFromh)(x1From, x2From, x3From);

	(*this->localDistributionsToh2)(D3Q27System::ET_E, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_E, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_N, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_N, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_T, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_T, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_NE, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_NE, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TE, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TE, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TN, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TN, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TNE, x1To, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TNE, x1From, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
    (*this->localDistributionsToh2)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) =
        (*this->localDistributionsFromh2)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_W, x1To + 1, x2To, x3To) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_S, x1To, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_B, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh2)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFromh2)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

    (*this->zeroDistributionsToh2)(x1To, x2To, x3To) = (*this->zeroDistributionsFromh2)(x1From, x2From, x3From);

    (*this->pressureTo)(x1To, x2To, x3To) = (*this->pressureFrom)(x1From, x2From, x3From);
}
#endif