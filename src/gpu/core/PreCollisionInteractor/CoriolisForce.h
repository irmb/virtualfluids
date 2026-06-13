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
//! \addtogroup gpu_PreCollisionInteractor PreCollisionInteractor
//! \{
//! \author Henry Korb
//! \date 27/02/2024
//=======================================================================================

#ifndef CORIOLIS_FORCE_H_
#define CORIOLIS_FORCE_H_

#include "Parameter/Parameter.h"
#include "PreCollisionInteractor.h"

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <logger/Logger.h>
#include <stdexcept>

namespace vf::gpu {

//!\brief Interactor to compute Coriolis force. All parameters in SI units.
class CoriolisForce : public PreCollisionInteractor
{
public:
    CoriolisForce(const SPtr<Parameter>& parameter, SPtr<CudaMemoryManager> cudaMemoryManager, real geostrophicWindX,
                  real geostrophicWindY, real coriolisParameter)
        : geostrophicWindX(geostrophicWindX),
          geostrophicWindY(geostrophicWindY),
          coriolisParameter(coriolisParameter),
          PreCollisionInteractor(parameter, std::move(cudaMemoryManager))
    {
        VF_LOG_INFO("using Coriolis Force with geostrophic wind vector ({},{}) m/s and coriolis parameter {} 1/s", geostrophicWindX,
                    geostrophicWindY, coriolisParameter);
        if (!para->getIsBodyForce())
            throw std::runtime_error("Coriolis force needs body force.");
        para->setAllNodesAllFeatures(true);
    }

    void init() override {};
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* /**/) override {};
    ~CoriolisForce() override = default;

private:
    const real geostrophicWindX, geostrophicWindY, coriolisParameter;
};

}

#endif //! \}