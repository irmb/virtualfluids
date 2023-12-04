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

#include "FullDirectConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "EsoSplit.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

//! \brief   Exchange data between blocks. 
//! \details Connector send and receive full distributions between two blocks in shared memory.

class TwoDistributionsFullDirectConnector : public FullDirectConnector
{
public:
    TwoDistributionsFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);
    void init() override;

    void exchangeData() override
    {
        FullDirectConnector::exchangeData();
    }

protected:
    inline void updatePointers() override;
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

    SPtr<EsoTwist3D>  fFrom, hFrom;
    SPtr<EsoTwist3D>  fTo, hTo;
};
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullDirectConnector::updatePointers()
{
    localDistributionsFromf = dynamicPointerCast<EsoSplit>(this->fFrom)->getLocalDistributions();
    nonLocalDistributionsFromf = dynamicPointerCast<EsoSplit>(this->fFrom)->getNonLocalDistributions();
    zeroDistributionsFromf = dynamicPointerCast<EsoSplit>(this->fFrom)->getZeroDistributions();

    localDistributionsTof    = dynamicPointerCast<EsoSplit>(this->fTo)->getLocalDistributions();
    nonLocalDistributionsTof = dynamicPointerCast<EsoSplit>(this->fTo)->getNonLocalDistributions();
    zeroDistributionsTof     = dynamicPointerCast<EsoSplit>(this->fTo)->getZeroDistributions();

    localDistributionsFromh = dynamicPointerCast<EsoSplit>(this->hFrom)->getLocalDistributions();
    nonLocalDistributionsFromh = dynamicPointerCast<EsoSplit>(this->hFrom)->getNonLocalDistributions();
    zeroDistributionsFromh = dynamicPointerCast<EsoSplit>(this->hFrom)->getZeroDistributions();

    localDistributionsToh    = dynamicPointerCast<EsoSplit>(this->hTo)->getLocalDistributions();
    nonLocalDistributionsToh = dynamicPointerCast<EsoSplit>(this->hTo)->getNonLocalDistributions();
    zeroDistributionsToh     = dynamicPointerCast<EsoSplit>(this->hTo)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{
    using namespace vf::lbm::dir;

    (*this->localDistributionsTof)(eP00, x1To, x2To, x3To) = (*this->localDistributionsFromf)(eP00, x1From, x2From, x3From);
    (*this->localDistributionsTof)(e0P0, x1To, x2To, x3To) = (*this->localDistributionsFromf)(e0P0, x1From, x2From, x3From);
    (*this->localDistributionsTof)(e00P, x1To, x2To, x3To) = (*this->localDistributionsFromf)(e00P, x1From, x2From, x3From);
    (*this->localDistributionsTof)(ePP0, x1To, x2To, x3To) = (*this->localDistributionsFromf)(ePP0, x1From, x2From, x3From);
    (*this->localDistributionsTof)(eMP0, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(eMP0, x1From + 1, x2From, x3From);
    (*this->localDistributionsTof)(eP0P, x1To, x2To, x3To) = (*this->localDistributionsFromf)(eP0P, x1From, x2From, x3From);
    (*this->localDistributionsTof)(eM0P, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(eM0P, x1From + 1, x2From, x3From);
    (*this->localDistributionsTof)(e0PP, x1To, x2To, x3To) = (*this->localDistributionsFromf)(e0PP, x1From, x2From, x3From);
    (*this->localDistributionsTof)(e0MP, x1To, x2To + 1, x3To) = (*this->localDistributionsFromf)(e0MP, x1From, x2From + 1, x3From);
    (*this->localDistributionsTof)(ePPP, x1To, x2To, x3To) = (*this->localDistributionsFromf)(ePPP, x1From, x2From, x3From);
    (*this->localDistributionsTof)(eMPP, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(eMPP, x1From + 1, x2From, x3From);
    (*this->localDistributionsTof)(ePMP, x1To, x2To + 1, x3To) = (*this->localDistributionsFromf)(ePMP, x1From, x2From + 1, x3From);
    (*this->localDistributionsTof)(eMMP, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFromf)(eMMP, x1From + 1, x2From + 1, x3From);

    (*this->nonLocalDistributionsTof)(eM00, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFromf)(eM00, x1From + 1, x2From, x3From);
    (*this->nonLocalDistributionsTof)(e0M0, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(e0M0, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsTof)(e00M, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(e00M, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTof)(eMM0, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(eMM0, x1From + 1, x2From + 1, x3From);
    (*this->nonLocalDistributionsTof)(ePM0, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(ePM0, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsTof)(eM0M, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(eM0M, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsTof)(eP0M, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(eP0M, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTof)(e0MM, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(e0MM, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTof)(e0PM, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(e0PM, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTof)(eMMM, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(eMMM, x1From + 1, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTof)(ePMM, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(ePMM, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTof)(eMPM, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(eMPM, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsTof)(ePPM, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(ePPM, x1From, x2From, x3From + 1);

    (*this->zeroDistributionsTof)(x1To, x2To, x3To) = (*this->zeroDistributionsFromf)(x1From, x2From, x3From);


    (*this->localDistributionsToh)(eP00, x1To, x2To, x3To) = (*this->localDistributionsFromh)(eP00, x1From, x2From, x3From);
    (*this->localDistributionsToh)(e0P0, x1To, x2To, x3To) = (*this->localDistributionsFromh)(e0P0, x1From, x2From, x3From);
    (*this->localDistributionsToh)(e00P, x1To, x2To, x3To) = (*this->localDistributionsFromh)(e00P, x1From, x2From, x3From);
    (*this->localDistributionsToh)(ePP0, x1To, x2To, x3To) = (*this->localDistributionsFromh)(ePP0, x1From, x2From, x3From);
    (*this->localDistributionsToh)(eMP0, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(eMP0, x1From + 1, x2From, x3From);
    (*this->localDistributionsToh)(eP0P, x1To, x2To, x3To) = (*this->localDistributionsFromh)(eP0P, x1From, x2From, x3From);
    (*this->localDistributionsToh)(eM0P, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(eM0P, x1From + 1, x2From, x3From);
    (*this->localDistributionsToh)(e0PP, x1To, x2To, x3To) = (*this->localDistributionsFromh)(e0PP, x1From, x2From, x3From);
    (*this->localDistributionsToh)(e0MP, x1To, x2To + 1, x3To) = (*this->localDistributionsFromh)(e0MP, x1From, x2From + 1, x3From);
    (*this->localDistributionsToh)(ePPP, x1To, x2To, x3To) = (*this->localDistributionsFromh)(ePPP, x1From, x2From, x3From);
    (*this->localDistributionsToh)(eMPP, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(eMPP, x1From + 1, x2From, x3From);
    (*this->localDistributionsToh)(ePMP, x1To, x2To + 1, x3To) = (*this->localDistributionsFromh)(ePMP, x1From, x2From + 1, x3From);
    (*this->localDistributionsToh)(eMMP, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFromh)(eMMP, x1From + 1, x2From + 1, x3From);

    (*this->nonLocalDistributionsToh)(eM00, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFromh)(eM00, x1From + 1, x2From, x3From);
    (*this->nonLocalDistributionsToh)(e0M0, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(e0M0, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsToh)(e00M, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(e00M, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh)(eMM0, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(eMM0, x1From + 1, x2From + 1, x3From);
    (*this->nonLocalDistributionsToh)(ePM0, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(ePM0, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsToh)(eM0M, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(eM0M, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh)(eP0M, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(eP0M, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh)(e0MM, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(e0MM, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsToh)(e0PM, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(e0PM, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh)(eMMM, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(eMMM, x1From + 1, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsToh)(ePMM, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(ePMM, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsToh)(eMPM, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(eMPM, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsToh)(ePPM, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(ePPM, x1From, x2From, x3From + 1);

    (*this->zeroDistributionsToh)(x1To, x2To, x3To) = (*this->zeroDistributionsFromh)(x1From, x2From, x3From);
}
#endif