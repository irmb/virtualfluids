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
//! \file LiggghtsCouplingLBMKernel.h
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef LiggghtsCouplingLBMKernel_h
#define LiggghtsCouplingLBMKernel_h

#include "LBMKernel.h"
#include "IBdynamicsParticleData.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class LiggghtsCouplingLBMKernel : public LBMKernel
{
public:
    virtual ~LiggghtsCouplingLBMKernel() = default;

    CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>::CbArray3DPtr getParticleData()
    {
        return particleData;
    };
    void setParticleData(CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>::CbArray3DPtr particleData)
    {
        this->particleData = particleData;
    };

    
 protected:
    //void collisionOperator();
    CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>::CbArray3DPtr particleData;

    //CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
    //CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
    //CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr restDistributionsF;
};

#endif