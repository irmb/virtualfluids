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
//! \ingroup gpu_core core
//! \{
#ifndef ActuatorFarmStandalone_H
#define ActuatorFarmStandalone_H

#include "ActuatorFarm.h"
#include "basics/DataTypes.h"

class ActuatorFarmStandalone : public ActuatorFarm
{
public:
    ActuatorFarmStandalone(
        const real diameter,
        const uint numberOfNodesPerBlade,
        const std::vector<real> turbinePositionsX,
        const std::vector<real> turbinePositionsY,
        const std::vector<real> turbinePositionsZ,
        const std::vector<real> rotorSpeeds,
        const real density,
        const real smearingWidth,
        const int level,
        const real deltaT,
        const real deltaX
    ) : rotorSpeeds(rotorSpeeds),
        ActuatorFarm(diameter, computeBladeRadii(diameter, numberOfNodesPerBlade), turbinePositionsX, turbinePositionsY, turbinePositionsZ, density, smearingWidth, level, deltaT, deltaX, true)
    {}

    ~ActuatorFarmStandalone() = default;

    void updateForcesAndCoordinates() override;
    static std::vector<real> computeBladeRadii(const real diameter, const uint numberOfNodesPerBlade);
    
private:
    std::vector<real> rotorSpeeds;

};

#endif

//! \}
