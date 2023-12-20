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

#ifndef OneDistributionFullDirectConnector_H
#define OneDistributionFullDirectConnector_H

#include "Block3D.h"
#include "D3Q27System.h"
#include "FullDirectConnector.h"
#include "EsoSplit.h"
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
    inline void updatePointers() override;

    void exchangeData() override
    {
        FullDirectConnector::exchangeData();
    }
    inline void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To) override;

private:
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFrom;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFrom;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsFrom;

    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTo;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTo;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsTo;

    SPtr<EsoTwist3D> fFrom;
    SPtr<EsoTwist3D> fTo;
};
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullDirectConnector::updatePointers()
{
    localDistributionsFrom = dynamicPointerCast<EsoSplit>(this->fFrom)->getLocalDistributions();
    nonLocalDistributionsFrom = dynamicPointerCast<EsoSplit>(this->fFrom)->getNonLocalDistributions();
    zeroDistributionsFrom = dynamicPointerCast<EsoSplit>(this->fFrom)->getZeroDistributions();

    localDistributionsTo = dynamicPointerCast<EsoSplit>(this->fTo)->getLocalDistributions();
    nonLocalDistributionsTo = dynamicPointerCast<EsoSplit>(this->fTo)->getNonLocalDistributions();
    zeroDistributionsTo     = dynamicPointerCast<EsoSplit>(this->fTo)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{
    using namespace vf::lbm::dir;

    (*this->localDistributionsTo)(eP00, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(eP00, x1From, x2From, x3From);
    (*this->localDistributionsTo)(e0P0, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(e0P0, x1From, x2From, x3From);
    (*this->localDistributionsTo)(e00P, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(e00P, x1From, x2From, x3From);
    (*this->localDistributionsTo)(ePP0, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(ePP0, x1From, x2From, x3From);
    (*this->localDistributionsTo)(eMP0, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFrom)(eMP0, x1From + 1, x2From, x3From);
    (*this->localDistributionsTo)(eP0P, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(eP0P, x1From, x2From, x3From);
    (*this->localDistributionsTo)(eM0P, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFrom)(eM0P, x1From + 1, x2From, x3From);
    (*this->localDistributionsTo)(e0PP, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(e0PP, x1From, x2From, x3From);
    (*this->localDistributionsTo)(e0MP, x1To, x2To + 1, x3To) =
        (*this->localDistributionsFrom)(e0MP, x1From, x2From + 1, x3From);
    (*this->localDistributionsTo)(ePPP, x1To, x2To, x3To) =
        (*this->localDistributionsFrom)(ePPP, x1From, x2From, x3From);
    (*this->localDistributionsTo)(eMPP, x1To + 1, x2To, x3To) =
        (*this->localDistributionsFrom)(eMPP, x1From + 1, x2From, x3From);
    (*this->localDistributionsTo)(ePMP, x1To, x2To + 1, x3To) =
        (*this->localDistributionsFrom)(ePMP, x1From, x2From + 1, x3From);
    (*this->localDistributionsTo)(eMMP, x1To + 1, x2To + 1, x3To) =
        (*this->localDistributionsFrom)(eMMP, x1From + 1, x2From + 1, x3From);

    (*this->nonLocalDistributionsTo)(eM00, x1To + 1, x2To, x3To) =
        (*this->nonLocalDistributionsFrom)(eM00, x1From + 1, x2From, x3From);
    (*this->nonLocalDistributionsTo)(e0M0, x1To, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFrom)(e0M0, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsTo)(e00M, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(e00M, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(eMM0, x1To + 1, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFrom)(eMM0, x1From + 1, x2From + 1, x3From);
    (*this->nonLocalDistributionsTo)(ePM0, x1To, x2To + 1, x3To) =
        (*this->nonLocalDistributionsFrom)(ePM0, x1From, x2From + 1, x3From);
    (*this->nonLocalDistributionsTo)(eM0M, x1To + 1, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(eM0M, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(eP0M, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(eP0M, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(e0MM, x1To, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(e0MM, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTo)(e0PM, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(e0PM, x1From, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(eMMM, x1To + 1, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(eMMM, x1From + 1, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTo)(ePMM, x1To, x2To + 1, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(ePMM, x1From, x2From + 1, x3From + 1);
    (*this->nonLocalDistributionsTo)(eMPM, x1To + 1, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(eMPM, x1From + 1, x2From, x3From + 1);
    (*this->nonLocalDistributionsTo)(ePPM, x1To, x2To, x3To + 1) =
        (*this->nonLocalDistributionsFrom)(ePPM, x1From, x2From, x3From + 1);

    (*this->zeroDistributionsTo)(x1To, x2To, x3To) = (*this->zeroDistributionsFrom)(x1From, x2From, x3From);
}
#endif

//! \}
