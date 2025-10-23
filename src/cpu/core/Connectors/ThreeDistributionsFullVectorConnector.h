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

#ifndef ThreeDistributionsFullVectorConnector_H
#define ThreeDistributionsFullVectorConnector_H

#include <vector>

#include "FullVectorConnector.h"
#include "D3Q27System.h"
#include "EsoSplit.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class EsoTwist3D;
class Block3D;

//! \brief   Exchange data between blocks.
//! \details Connector send and receive full distributions between two blocks in distributed memory.
class ThreeDistributionsFullVectorConnector : public FullVectorConnector
{
public:
    ThreeDistributionsFullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir);

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
   CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr splitA;
   CbArray4D <real, IndexerX4X3X2X1>::CbArray4DPtr splitB;
   CbArray3D <real, IndexerX3X2X1>::CbArray3DPtr   split0;

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
    splitA    = dynamicPointerCast<EsoSplit>(this->fDis)->getSplitA();
    splitB = dynamicPointerCast<EsoSplit>(this->fDis)->getSplitB();
    split0     = dynamicPointerCast<EsoSplit>(this->fDis)->getSplit0();

    localHdistributions    = dynamicPointerCast<EsoSplit>(this->hDis)->getSplitA();
    nonLocalHdistributions = dynamicPointerCast<EsoSplit>(this->hDis)->getSplitB();
    zeroHdistributions     = dynamicPointerCast<EsoSplit>(this->hDis)->getSplit0();

    localH2distributions    = dynamicPointerCast<EsoSplit>(this->h2Dis)->getSplitA();
    nonLocalH2distributions = dynamicPointerCast<EsoSplit>(this->h2Dis)->getSplitB();
    zeroH2distributions     = dynamicPointerCast<EsoSplit>(this->h2Dis)->getSplit0();
}
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsFullVectorConnector::fillData(vector_type& sdata, int& index, int x1, int x2, int x3)
{
   using namespace vf::lbm::dir;

   sdata[index++] = (*this->splitA)(eP00, x1, x2, x3);
   sdata[index++] = (*this->splitA)(e0P0, x1, x2, x3);
   sdata[index++] = (*this->splitA)(e00P, x1, x2, x3);
   sdata[index++] = (*this->splitA)(ePP0, x1, x2, x3);
   sdata[index++] = (*this->splitA)(eMP0, x1 + 1, x2, x3);
   sdata[index++] = (*this->splitA)(eP0P, x1, x2, x3);
   sdata[index++] = (*this->splitA)(eM0P, x1 + 1, x2, x3);
   sdata[index++] = (*this->splitA)(e0PP, x1, x2, x3);
   sdata[index++] = (*this->splitA)(e0MP, x1, x2 + 1, x3);
   sdata[index++] = (*this->splitA)(ePPP, x1, x2, x3);
   sdata[index++] = (*this->splitA)(eMPP, x1 + 1, x2, x3);
   sdata[index++] = (*this->splitA)(ePMP, x1, x2 + 1, x3);
   sdata[index++] = (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->splitB)(eM00, x1 + 1, x2, x3);
   sdata[index++] = (*this->splitB)(e0M0, x1, x2 + 1, x3);
   sdata[index++] = (*this->splitB)(e00M, x1, x2, x3 + 1);
   sdata[index++] = (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->splitB)(ePM0, x1, x2 + 1, x3);
   sdata[index++] = (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->splitB)(eP0M, x1, x2, x3 + 1);
   sdata[index++] = (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->splitB)(e0PM, x1, x2, x3 + 1);
   sdata[index++] = (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->splitB)(ePPM, x1, x2, x3 + 1);

   sdata[index++] = (*this->split0)(x1, x2, x3);


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

   sdata[index++] = (*this->localH2distributions)(eP00, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(e0P0, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(e00P, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(ePP0, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(eMP0, x1 + 1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(eP0P, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(eM0P, x1 + 1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(e0PP, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(e0MP, x1, x2 + 1, x3);
   sdata[index++] = (*this->localH2distributions)(ePPP, x1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(eMPP, x1 + 1, x2, x3);
   sdata[index++] = (*this->localH2distributions)(ePMP, x1, x2 + 1, x3);
   sdata[index++] = (*this->localH2distributions)(eMMP, x1 + 1, x2 + 1, x3);

   sdata[index++] = (*this->nonLocalH2distributions)(eM00, x1 + 1, x2, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(e0M0, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(e00M, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(eMM0, x1 + 1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(ePM0, x1, x2 + 1, x3);
   sdata[index++] = (*this->nonLocalH2distributions)(eM0M, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(eP0M, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(e0MM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(e0PM, x1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(ePMM, x1, x2 + 1, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(eMPM, x1 + 1, x2, x3 + 1);
   sdata[index++] = (*this->nonLocalH2distributions)(ePPM, x1, x2, x3 + 1);

   sdata[index++] = (*this->zeroH2distributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
inline void ThreeDistributionsFullVectorConnector::distributeData(vector_type& rdata, int& index, int x1, int x2, int x3)
{
   using namespace vf::lbm::dir;

   (*this->splitA)(eP00, x1, x2, x3) = rdata[index++];
   (*this->splitA)(e0P0, x1, x2, x3) = rdata[index++];
   (*this->splitA)(e00P, x1, x2, x3) = rdata[index++];
   (*this->splitA)(ePP0, x1, x2, x3) = rdata[index++];
   (*this->splitA)(eMP0, x1 + 1, x2, x3) = rdata[index++];
   (*this->splitA)(eP0P, x1, x2, x3) = rdata[index++];
   (*this->splitA)(eM0P, x1 + 1, x2, x3) = rdata[index++];
   (*this->splitA)(e0PP, x1, x2, x3) = rdata[index++];
   (*this->splitA)(e0MP, x1, x2 + 1, x3) = rdata[index++];
   (*this->splitA)(ePPP, x1, x2, x3) = rdata[index++];
   (*this->splitA)(eMPP, x1 + 1, x2, x3) = rdata[index++];
   (*this->splitA)(ePMP, x1, x2 + 1, x3) = rdata[index++];
   (*this->splitA)(eMMP, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->splitB)(eM00, x1 + 1, x2, x3) = rdata[index++];
   (*this->splitB)(e0M0, x1, x2 + 1, x3) = rdata[index++];
   (*this->splitB)(e00M, x1, x2, x3 + 1) = rdata[index++];
   (*this->splitB)(eMM0, x1 + 1, x2 + 1, x3) = rdata[index++];
   (*this->splitB)(ePM0, x1, x2 + 1, x3) = rdata[index++];
   (*this->splitB)(eM0M, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->splitB)(eP0M, x1, x2, x3 + 1) = rdata[index++];
   (*this->splitB)(e0MM, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->splitB)(e0PM, x1, x2, x3 + 1) = rdata[index++];
   (*this->splitB)(eMMM, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->splitB)(ePMM, x1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->splitB)(eMPM, x1 + 1, x2, x3 + 1) = rdata[index++];
   (*this->splitB)(ePPM, x1, x2, x3 + 1) = rdata[index++];

   (*this->split0)(x1, x2, x3) = rdata[index++];

   
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

   (*this->localH2distributions)(eP00, x1, x2, x3)           = rdata[index++];
   (*this->localH2distributions)(e0P0, x1, x2, x3)           = rdata[index++];
   (*this->localH2distributions)(e00P, x1, x2, x3)           = rdata[index++];
   (*this->localH2distributions)(ePP0, x1, x2, x3)          = rdata[index++];
   (*this->localH2distributions)(eMP0, x1 + 1, x2, x3)      = rdata[index++];
   (*this->localH2distributions)(eP0P, x1, x2, x3)          = rdata[index++];
   (*this->localH2distributions)(eM0P, x1 + 1, x2, x3)      = rdata[index++];
   (*this->localH2distributions)(e0PP, x1, x2, x3)          = rdata[index++];
   (*this->localH2distributions)(e0MP, x1, x2 + 1, x3)      = rdata[index++];
   (*this->localH2distributions)(ePPP, x1, x2, x3)         = rdata[index++];
   (*this->localH2distributions)(eMPP, x1 + 1, x2, x3)     = rdata[index++];
   (*this->localH2distributions)(ePMP, x1, x2 + 1, x3)     = rdata[index++];
   (*this->localH2distributions)(eMMP, x1 + 1, x2 + 1, x3) = rdata[index++];

   (*this->nonLocalH2distributions)(eM00, x1 + 1, x2, x3)           = rdata[index++];
   (*this->nonLocalH2distributions)(e0M0, x1, x2 + 1, x3)           = rdata[index++];
   (*this->nonLocalH2distributions)(e00M, x1, x2, x3 + 1)           = rdata[index++];
   (*this->nonLocalH2distributions)(eMM0, x1 + 1, x2 + 1, x3)      = rdata[index++];
   (*this->nonLocalH2distributions)(ePM0, x1, x2 + 1, x3)          = rdata[index++];
   (*this->nonLocalH2distributions)(eM0M, x1 + 1, x2, x3 + 1)      = rdata[index++];
   (*this->nonLocalH2distributions)(eP0M, x1, x2, x3 + 1)          = rdata[index++];
   (*this->nonLocalH2distributions)(e0MM, x1, x2 + 1, x3 + 1)      = rdata[index++];
   (*this->nonLocalH2distributions)(e0PM, x1, x2, x3 + 1)          = rdata[index++];
   (*this->nonLocalH2distributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
   (*this->nonLocalH2distributions)(ePMM, x1, x2 + 1, x3 + 1)     = rdata[index++];
   (*this->nonLocalH2distributions)(eMPM, x1 + 1, x2, x3 + 1)     = rdata[index++];
   (*this->nonLocalH2distributions)(ePPM, x1, x2, x3 + 1)         = rdata[index++];

   (*this->zeroH2distributions)(x1, x2, x3) = rdata[index++];
}


#endif 


//! \}
