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
//! \file TwoDistributionsFullVectorConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef TwoDistributionsFullVectorConnector_H
#define TwoDistributionsFullVectorConnector_H

#include <vector>

#include "RemoteBlock3DConnector.h"
#include "D3Q27System.h"
//#include "Block3D.h"
//#include "LBMKernel.h"
//#include "EsoTwist3D.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class EsoTwist3D;
class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
class TwoDistributionsFullVectorConnector : public RemoteBlock3DConnector
{
public:
   TwoDistributionsFullVectorConnector(SPtr<Block3D> block
      , VectorTransmitterPtr sender
      , VectorTransmitterPtr receiver
      , int sendDir);

   void init();

   void fillSendVectors();
   void distributeReceiveVectors();

protected:
   inline void fillData(vector_type& sdata, int& index, int x1, int x2, int x3);
   inline void distributeData(vector_type& rdata, int& index, int x1, int x2, int x3);
private:
   int maxX1;
   int maxX2;
   int maxX3;

   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   SPtr<EsoTwist3D>  fDis;

   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localHdistributions;
   CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalHdistributions;
   CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroHdistributions;


   SPtr<EsoTwist3D>  hDis;

};

//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullVectorConnector::fillData(vector_type& sdata, int& index, int x1, int x2, int x3)
{
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
   sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroDistributions)(x1, x2, x3);


   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_E, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_N, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_T, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_NE, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TE, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TN, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TNE, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
   sdata[index++] = (*this->localHdistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroHdistributions)(x1, x2, x3);

}
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullVectorConnector::distributeData(vector_type& rdata, int& index, int x1, int x2, int x3)
{
   (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = rdata[index++];
   (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = rdata[index++];

   (*this->zeroDistributions)(x1, x2, x3) = rdata[index++];

   
   (*this->localHdistributions)(D3Q27System::ET_E, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_N, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_T, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_NE, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TE, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TN, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TNE, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3) = rdata[index++];
   (*this->localHdistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalHdistributions)(D3Q27System::ET_W, x1 + 1, x2, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_S, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_B, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1) = rdata[index++];

   (*this->zeroHdistributions)(x1, x2, x3) = rdata[index++];
}


#endif //D3Q27VECTORCONNECTOR_H

