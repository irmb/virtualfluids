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
//! \file TwoDistributionsFullDirectConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef TwoDistributionsFullDirectConnector_H
#define TwoDistributionsFullDirectConnector_H

#include "LocalBlock3DConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"
#include "BCArray3D.h"

class EsoTwist3D;

//! \brief   Exchange data between blocks. 
//! \details Connector send and receive full distributions between two blocks in shared memory.

class TwoDistributionsFullDirectConnector : public LocalBlock3DConnector
{
public:
	TwoDistributionsFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);
	void init();
	void sendVectors();
	//void sendCellVectors();

protected:
	inline void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);
private:
	int maxX1;
	int maxX2;
	int maxX3;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromf;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromf;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFromf;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTof;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTof;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsTof;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromh;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromh;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFromh;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsToh;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsToh;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsToh;

	SPtr<EsoTwist3D>  fFrom, hFrom;
	SPtr<EsoTwist3D>  fTo, hTo;
};
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
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
}
#endif