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

#include "FullVectorConnector.h"
#include "D3Q27System.h"
#include "EsoSplit.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class EsoTwist3D;
class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
class TwoDistributionsFullVectorConnector : public FullVectorConnector
{
public:
    TwoDistributionsFullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir);

    void init() override;

    void fillData() override
    {
        FullVectorConnector::fillData();
    }

    void distributeData() override
    {
        FullVectorConnector::distributeData();
    }

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

};
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullVectorConnector::updatePointers()
{
    localDistributions    = dynamicPointerCast<EsoSplit>(this->fDis)->getLocalDistributions();
    nonLocalDistributions = dynamicPointerCast<EsoSplit>(this->fDis)->getNonLocalDistributions();
    zeroDistributions     = dynamicPointerCast<EsoSplit>(this->fDis)->getZeroDistributions();

    localHdistributions    = dynamicPointerCast<EsoSplit>(this->hDis)->getLocalDistributions();
    nonLocalHdistributions = dynamicPointerCast<EsoSplit>(this->hDis)->getNonLocalDistributions();
    zeroHdistributions     = dynamicPointerCast<EsoSplit>(this->hDis)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullVectorConnector::fillData(vector_type& sdata, int& index, int x1, int x2, int x3)
{
   using namespace vf::lbm::dir;

   sdata[index++] = (*this->localDistributions)(eP00, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(e0P0, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(e00P, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(ePP0, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(eMP0, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(eP0P, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(eM0P, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(e0PP, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(e0MP, x1, x2 + 1, x3);
   sdata[index++] = (*this->localDistributions)(ePPP, x1, x2, x3);
   sdata[index++] = (*this->localDistributions)(eMPP, x1 + 1, x2, x3);
   sdata[index++] = (*this->localDistributions)(ePMP, x1, x2 + 1, x3);
   sdata[index++] = (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroDistributions)(x1, x2, x3);


   sdata[index++] = (*this->localHdistributions)(eP00, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(e0P0, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(e00P, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(ePP0, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(eMP0, x1 + 1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(eP0P, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(eM0P, x1 + 1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(e0PP, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(e0MP, x1, x2 + 1, x3);
   sdata[index++] = (*this->localHdistributions)(ePPP, x1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(eMPP, x1 + 1, x2, x3);
   sdata[index++] = (*this->localHdistributions)(ePMP, x1, x2 + 1, x3);
   sdata[index++] = (*this->localHdistributions)(eMMP, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalHdistributions)(eM00, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(e0M0, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(e00M, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(eMM0, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(ePM0, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalHdistributions)(eM0M, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(eP0M, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(e0MM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(e0PM, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(ePMM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(eMPM, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalHdistributions)(ePPM, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroHdistributions)(x1, x2, x3);

}
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullVectorConnector::distributeData(vector_type& rdata, int& index, int x1, int x2, int x3)
{
   using namespace vf::lbm::dir;

   (*this->localDistributions)(eP00, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(e0P0, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(e00P, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(ePP0, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(eMP0, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(eP0P, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(eM0P, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(e0PP, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(e0MP, x1, x2 + 1, x3) = rdata[index++];
   (*this->localDistributions)(ePPP, x1, x2, x3) = rdata[index++];
   (*this->localDistributions)(eMPP, x1 + 1, x2, x3) = rdata[index++];
   (*this->localDistributions)(ePMP, x1, x2 + 1, x3) = rdata[index++];
   (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3) = rdata[index++];
   (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1) = rdata[index++];

   (*this->zeroDistributions)(x1, x2, x3) = rdata[index++];

   
   (*this->localHdistributions)(eP00, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(e0P0, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(e00P, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(ePP0, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(eMP0, x1 + 1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(eP0P, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(eM0P, x1 + 1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(e0PP, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(e0MP, x1, x2 + 1, x3) = rdata[index++];
   (*this->localHdistributions)(ePPP, x1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(eMPP, x1 + 1, x2, x3) = rdata[index++];
   (*this->localHdistributions)(ePMP, x1, x2 + 1, x3) = rdata[index++];
   (*this->localHdistributions)(eMMP, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalHdistributions)(eM00, x1 + 1, x2, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(e0M0, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(e00M, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(eMM0, x1 + 1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(ePM0, x1, x2 + 1, x3) = rdata[index++];
   (*this->nonLocalHdistributions)(eM0M, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(eP0M, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(e0MM, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(e0PM, x1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(ePMM, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(eMPM, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->nonLocalHdistributions)(ePPM, x1, x2, x3 + 1) = rdata[index++];

   (*this->zeroHdistributions)(x1, x2, x3) = rdata[index++];
}


#endif 

