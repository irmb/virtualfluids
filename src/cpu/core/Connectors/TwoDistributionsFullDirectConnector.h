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
//! \addtogroup cpu_Connectors Connectors
//! \ingroup cpu_core core
//! \{
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
    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr fSplitAFrom;
    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr fSplitBFrom;
    CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   fSplit0from;

    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr fSplitATo;
    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr fSplitBTo;
    CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   fSplit0To;

    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr hSplitAFrom;
    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr hSplitBFrom;
    CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   hSplit0From;

    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr hSplitATo;
    CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr hSplitBTo;
    CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   hSplit0To;

    SPtr<EsoTwist3D>  fFrom, hFrom;
    SPtr<EsoTwist3D>  fTo, hTo;
};
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullDirectConnector::updatePointers()
{
    fSplitAFrom= dynamicPointerCast<EsoSplit>(this->fFrom)->getSplitA();
    fSplitBFrom = dynamicPointerCast<EsoSplit>(this->fFrom)->getSplitB();
    fSplit0from = dynamicPointerCast<EsoSplit>(this->fFrom)->getSplit0();

    fSplitATo   = dynamicPointerCast<EsoSplit>(this->fTo)->getSplitA();
    fSplitBTo= dynamicPointerCast<EsoSplit>(this->fTo)->getSplitB();
    fSplit0To     = dynamicPointerCast<EsoSplit>(this->fTo)->getSplit0();

    hSplitAFrom = dynamicPointerCast<EsoSplit>(this->hFrom)->getSplitA();
    hSplitBFrom = dynamicPointerCast<EsoSplit>(this->hFrom)->getSplitB();
    hSplit0From = dynamicPointerCast<EsoSplit>(this->hFrom)->getSplit0();

    hSplitATo    = dynamicPointerCast<EsoSplit>(this->hTo)->getSplitA();
    hSplitBTo = dynamicPointerCast<EsoSplit>(this->hTo)->getSplitB();
    hSplit0To     = dynamicPointerCast<EsoSplit>(this->hTo)->getSplit0();
}
//////////////////////////////////////////////////////////////////////////
inline void TwoDistributionsFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{
    using namespace vf::lbm::dir;

    (*this->fSplitATo)(eP00, x1To, x2To, x3To) = (*this->fSplitAFrom)(eP00, x1From, x2From, x3From);
    (*this->fSplitATo)(e0P0, x1To, x2To, x3To) = (*this->fSplitAFrom)(e0P0, x1From, x2From, x3From);
    (*this->fSplitATo)(e00P, x1To, x2To, x3To) = (*this->fSplitAFrom)(e00P, x1From, x2From, x3From);
    (*this->fSplitATo)(ePP0, x1To, x2To, x3To) = (*this->fSplitAFrom)(ePP0, x1From, x2From, x3From);
    (*this->fSplitATo)(eMP0, x1To + 1, x2To, x3To) = (*this->fSplitAFrom)(eMP0, x1From + 1, x2From, x3From);
    (*this->fSplitATo)(eP0P, x1To, x2To, x3To) = (*this->fSplitAFrom)(eP0P, x1From, x2From, x3From);
    (*this->fSplitATo)(eM0P, x1To + 1, x2To, x3To) = (*this->fSplitAFrom)(eM0P, x1From + 1, x2From, x3From);
    (*this->fSplitATo)(e0PP, x1To, x2To, x3To) = (*this->fSplitAFrom)(e0PP, x1From, x2From, x3From);
    (*this->fSplitATo)(e0MP, x1To, x2To + 1, x3To) = (*this->fSplitAFrom)(e0MP, x1From, x2From + 1, x3From);
    (*this->fSplitATo)(ePPP, x1To, x2To, x3To) = (*this->fSplitAFrom)(ePPP, x1From, x2From, x3From);
    (*this->fSplitATo)(eMPP, x1To + 1, x2To, x3To) = (*this->fSplitAFrom)(eMPP, x1From + 1, x2From, x3From);
    (*this->fSplitATo)(ePMP, x1To, x2To + 1, x3To) = (*this->fSplitAFrom)(ePMP, x1From, x2From + 1, x3From);
    (*this->fSplitATo)(eMMP, x1To + 1, x2To + 1, x3To) = (*this->fSplitAFrom)(eMMP, x1From + 1, x2From + 1, x3From);

    (*this->fSplitBTo)(eM00, x1To + 1, x2To, x3To) = (*this->fSplitBFrom)(eM00, x1From + 1, x2From, x3From);
    (*this->fSplitBTo)(e0M0, x1To, x2To + 1, x3To) = (*this->fSplitBFrom)(e0M0, x1From, x2From + 1, x3From);
    (*this->fSplitBTo)(e00M, x1To, x2To, x3To + 1) = (*this->fSplitBFrom)(e00M, x1From, x2From, x3From + 1);
    (*this->fSplitBTo)(eMM0, x1To + 1, x2To + 1, x3To) = (*this->fSplitBFrom)(eMM0, x1From + 1, x2From + 1, x3From);
    (*this->fSplitBTo)(ePM0, x1To, x2To + 1, x3To) = (*this->fSplitBFrom)(ePM0, x1From, x2From + 1, x3From);
    (*this->fSplitBTo)(eM0M, x1To + 1, x2To, x3To + 1) = (*this->fSplitBFrom)(eM0M, x1From + 1, x2From, x3From + 1);
    (*this->fSplitBTo)(eP0M, x1To, x2To, x3To + 1) = (*this->fSplitBFrom)(eP0M, x1From, x2From, x3From + 1);
    (*this->fSplitBTo)(e0MM, x1To, x2To + 1, x3To + 1) = (*this->fSplitBFrom)(e0MM, x1From, x2From + 1, x3From + 1);
    (*this->fSplitBTo)(e0PM, x1To, x2To, x3To + 1) = (*this->fSplitBFrom)(e0PM, x1From, x2From, x3From + 1);
    (*this->fSplitBTo)(eMMM, x1To + 1, x2To + 1, x3To + 1) = (*this->fSplitBFrom)(eMMM, x1From + 1, x2From + 1, x3From + 1);
    (*this->fSplitBTo)(ePMM, x1To, x2To + 1, x3To + 1) = (*this->fSplitBFrom)(ePMM, x1From, x2From + 1, x3From + 1);
    (*this->fSplitBTo)(eMPM, x1To + 1, x2To, x3To + 1) = (*this->fSplitBFrom)(eMPM, x1From + 1, x2From, x3From + 1);
    (*this->fSplitBTo)(ePPM, x1To, x2To, x3To + 1) = (*this->fSplitBFrom)(ePPM, x1From, x2From, x3From + 1);

    (*this->fSplit0To)(x1To, x2To, x3To) = (*this->fSplit0from)(x1From, x2From, x3From);


    (*this->hSplitATo)(eP00, x1To, x2To, x3To) = (*this->hSplitAFrom)(eP00, x1From, x2From, x3From);
    (*this->hSplitATo)(e0P0, x1To, x2To, x3To) = (*this->hSplitAFrom)(e0P0, x1From, x2From, x3From);
    (*this->hSplitATo)(e00P, x1To, x2To, x3To) = (*this->hSplitAFrom)(e00P, x1From, x2From, x3From);
    (*this->hSplitATo)(ePP0, x1To, x2To, x3To) = (*this->hSplitAFrom)(ePP0, x1From, x2From, x3From);
    (*this->hSplitATo)(eMP0, x1To + 1, x2To, x3To) = (*this->hSplitAFrom)(eMP0, x1From + 1, x2From, x3From);
    (*this->hSplitATo)(eP0P, x1To, x2To, x3To) = (*this->hSplitAFrom)(eP0P, x1From, x2From, x3From);
    (*this->hSplitATo)(eM0P, x1To + 1, x2To, x3To) = (*this->hSplitAFrom)(eM0P, x1From + 1, x2From, x3From);
    (*this->hSplitATo)(e0PP, x1To, x2To, x3To) = (*this->hSplitAFrom)(e0PP, x1From, x2From, x3From);
    (*this->hSplitATo)(e0MP, x1To, x2To + 1, x3To) = (*this->hSplitAFrom)(e0MP, x1From, x2From + 1, x3From);
    (*this->hSplitATo)(ePPP, x1To, x2To, x3To) = (*this->hSplitAFrom)(ePPP, x1From, x2From, x3From);
    (*this->hSplitATo)(eMPP, x1To + 1, x2To, x3To) = (*this->hSplitAFrom)(eMPP, x1From + 1, x2From, x3From);
    (*this->hSplitATo)(ePMP, x1To, x2To + 1, x3To) = (*this->hSplitAFrom)(ePMP, x1From, x2From + 1, x3From);
    (*this->hSplitATo)(eMMP, x1To + 1, x2To + 1, x3To) = (*this->hSplitAFrom)(eMMP, x1From + 1, x2From + 1, x3From);

    (*this->hSplitBTo)(eM00, x1To + 1, x2To, x3To) = (*this->hSplitBFrom)(eM00, x1From + 1, x2From, x3From);
    (*this->hSplitBTo)(e0M0, x1To, x2To + 1, x3To) = (*this->hSplitBFrom)(e0M0, x1From, x2From + 1, x3From);
    (*this->hSplitBTo)(e00M, x1To, x2To, x3To + 1) = (*this->hSplitBFrom)(e00M, x1From, x2From, x3From + 1);
    (*this->hSplitBTo)(eMM0, x1To + 1, x2To + 1, x3To) = (*this->hSplitBFrom)(eMM0, x1From + 1, x2From + 1, x3From);
    (*this->hSplitBTo)(ePM0, x1To, x2To + 1, x3To) = (*this->hSplitBFrom)(ePM0, x1From, x2From + 1, x3From);
    (*this->hSplitBTo)(eM0M, x1To + 1, x2To, x3To + 1) = (*this->hSplitBFrom)(eM0M, x1From + 1, x2From, x3From + 1);
    (*this->hSplitBTo)(eP0M, x1To, x2To, x3To + 1) = (*this->hSplitBFrom)(eP0M, x1From, x2From, x3From + 1);
    (*this->hSplitBTo)(e0MM, x1To, x2To + 1, x3To + 1) = (*this->hSplitBFrom)(e0MM, x1From, x2From + 1, x3From + 1);
    (*this->hSplitBTo)(e0PM, x1To, x2To, x3To + 1) = (*this->hSplitBFrom)(e0PM, x1From, x2From, x3From + 1);
    (*this->hSplitBTo)(eMMM, x1To + 1, x2To + 1, x3To + 1) = (*this->hSplitBFrom)(eMMM, x1From + 1, x2From + 1, x3From + 1);
    (*this->hSplitBTo)(ePMM, x1To, x2To + 1, x3To + 1) = (*this->hSplitBFrom)(ePMM, x1From, x2From + 1, x3From + 1);
    (*this->hSplitBTo)(eMPM, x1To + 1, x2To, x3To + 1) = (*this->hSplitBFrom)(eMPM, x1From + 1, x2From, x3From + 1);
    (*this->hSplitBTo)(ePPM, x1To, x2To, x3To + 1) = (*this->hSplitBFrom)(ePPM, x1From, x2From, x3From + 1);

    (*this->hSplit0To)(x1To, x2To, x3To) = (*this->hSplit0From)(x1From, x2From, x3From);
}
#endif
//! \}
