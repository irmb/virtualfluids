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
//! \file ThreeDistributionsFullVectorConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef ThreeDistributionsFullVectorConnector_H
#define ThreeDistributionsFullVectorConnector_H

#include <vector>

#include "FullVectorConnector.h"
#include "D3Q27System.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class EsoTwist3D;
class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
class ThreeDistributionsFullVectorConnector : public FullVectorConnector
{
public:
   ThreeDistributionsFullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir);

   void init() override;

protected:
   inline void updatePointers() override;
   inline void fillData(vector_type &sdata, int &index, int x1, int x2, int x3) override;
   inline void distributeData(vector_type &rdata, int &index, int x1, int x2, int x3) override;

private:
   CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   SPtr<EsoTwist3D>  fDis;

   CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr localHdistributions;
   CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalHdistributions;
   CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   zeroHdistributions;

   SPtr<EsoTwist3D>  hDis;

   CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localH2distributions;
   CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalH2distributions;
   CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroH2distributions;

   SPtr<EsoTwist3D> h2Dis;

};
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsFullVectorConnector::updatePointers()
{
    localDistributions    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getLocalDistributions();
    nonLocalDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getNonLocalDistributions();
    zeroDistributions     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getZeroDistributions();

    localHdistributions    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hDis)->getLocalDistributions();
    nonLocalHdistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hDis)->getNonLocalDistributions();
    zeroHdistributions     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hDis)->getZeroDistributions();

    localH2distributions    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->h2Dis)->getLocalDistributions();
    nonLocalH2distributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->h2Dis)->getNonLocalDistributions();
    zeroH2distributions     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->h2Dis)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsFullVectorConnector::fillData(vector_type& sdata, int& index, int x1, int x2, int x3)
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

   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_E, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_N, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_T, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_NE, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TE, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TN, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TNE, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
   sdata[index++] = (*this->localH2distributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroH2distributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsFullVectorConnector::distributeData(vector_type& rdata, int& index, int x1, int x2, int x3)
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

   (*this->localH2distributions)(D3Q27System::ET_E, x1, x2, x3)           = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_N, x1, x2, x3)           = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_T, x1, x2, x3)           = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_NE, x1, x2, x3)          = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_NW, x1 + 1, x2, x3)      = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TE, x1, x2, x3)          = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TW, x1 + 1, x2, x3)      = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TN, x1, x2, x3)          = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TS, x1, x2 + 1, x3)      = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TNE, x1, x2, x3)         = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3)     = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3)     = rdata[index++];
   (*this->localH2distributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalH2distributions)(D3Q27System::ET_W, x1 + 1, x2, x3)           = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_S, x1, x2 + 1, x3)           = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_B, x1, x2, x3 + 1)           = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3)      = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_SE, x1, x2 + 1, x3)          = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1)      = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BE, x1, x2, x3 + 1)          = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1)      = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BN, x1, x2, x3 + 1)          = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1)     = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1)     = rdata[index++];
   (*this->nonLocalH2distributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1)         = rdata[index++];

   (*this->zeroH2distributions)(x1, x2, x3) = rdata[index++];
}


#endif 

