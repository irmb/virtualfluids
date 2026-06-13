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

namespace vf::gpu {

class ActuatorFarmStandalone : public ActuatorFarm
{
public:
    ActuatorFarmStandalone(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, const real diameter,
                           const uint numberOfPointsPerBlade, const std::vector<real>& turbinePositionsX,
                           const std::vector<real>& turbinePositionsY, const std::vector<real>& turbinePositionsZ,
                           const std::vector<real>& rotorSpeeds, const real smearingWidth, const int level,
                           const std::optional<HubConfig>& hubConfig = std::nullopt,
                           const std::optional<TowerConfig>& towerConfig = std::nullopt,
                           const std::optional<real> hubDragCoeff = std::nullopt,
                           const std::optional<real> hubSkinFrictionCoeff = std::nullopt,
                           const std::optional<real> towerDragCoeff = std::nullopt,
                           const uint numberOfBlades = 3,
                           const std::optional<std::vector<real>>& bladeNormalCoefficients = std::nullopt)
        : rotorSpeeds(rotorSpeeds),
          hubDragCoeff(hubDragCoeff),
          hubSkinFrictionCoeff(hubSkinFrictionCoeff),
          towerDragCoeff(towerDragCoeff),
          bladeNormalCoefficients(bladeNormalCoefficients),
          ActuatorFarm(std::move(para), std::move(cudaMemoryManager), diameter,
                       computeBladeRadii(diameter, numberOfPointsPerBlade), turbinePositionsX, turbinePositionsY,
                       turbinePositionsZ, smearingWidth, level, true, numberOfBlades, hubConfig, towerConfig)
    {
        if (numberOfTurbines != rotorSpeeds.size())
            throw std::runtime_error("ActuatorFarmStandalone::ActuatorFarmStandalone: rotor speeds need to have same length "
                                     "as turbine positions!");
        if (bladeNormalCoefficients.has_value() && bladeNormalCoefficients->size() != numberOfPointsPerBlade)
            throw std::runtime_error("ActuatorFarmStandalone::ActuatorFarmStandalone: bladeNormalCoefficients need to have "
                                     "same length as numberOfPointsPerBlade.");
        
        if (numberOfTowerPointsPerTurbine > 0 && !towerDragCoeff.has_value())
            throw std::runtime_error(
                "Tower geometry is defined (tower points per turbine > 0), but towerDragCoeff is not provided!");

        if (numberOfHubPointsPerTurbine > 0)
        {
            if (!hubDragCoeff.has_value())
                throw std::runtime_error(
                    "Hub geometry is defined (hub points per turbine > 0), but hubDragCoeff is not provided!");

            if (!hubSkinFrictionCoeff.has_value())
                throw std::runtime_error("Hub geometry is defined (hub points per turbine > 0), but hubSkinFrictionCoeff is "
                                        "not provided!");
        }

        if (towerDragCoeff.has_value() && towerDragCoeff.value() > 0.0 && numberOfTowerPointsPerTurbine == 0)
            VF_LOG_WARNING("towerDragCoeff > 0 was provided, but no tower points per turbine exist → no tower forces will be applied.");

        if (hubDragCoeff.has_value() && hubDragCoeff.value() > 0.0 && numberOfHubPointsPerTurbine == 0)
            VF_LOG_WARNING("hubDragCoeff > 0 was provided, but no hub points per turbine exist → no hub forces will be applied.");

        if (hubSkinFrictionCoeff.has_value() && hubSkinFrictionCoeff.value() > 0.0 && numberOfHubPointsPerTurbine == 0)
            VF_LOG_WARNING("hubSkinFrictionCoeff > 0 was provided, but no hub points per turbine exist → no hub forces will be applied.");

        if (numberOfTowerPointsPerTurbine > 0 && towerDragCoeff.has_value() && towerDragCoeff.value() == 0.0)
            VF_LOG_WARNING("towerDragCoeff is set to 0 → no tower drag forces will be applied despite existing tower geometry.");

        if (numberOfHubPointsPerTurbine > 0)
        {
            if (hubDragCoeff.has_value() && hubDragCoeff.value() == 0.0)
                VF_LOG_WARNING("hubDragCoeff is set to 0 → no hub drag forces will be applied despite existing hub geometry.");

            if (hubSkinFrictionCoeff.has_value() && hubSkinFrictionCoeff.value() == 0.0)
                VF_LOG_WARNING("hubSkinFrictionCoeff is set to 0 → no hub skin friction forces will be applied despite existing hub geometry.");
        }

        VF_LOG_INFO("rotor speed [rad/s] = {}", this->rotorSpeeds[0]);

        if (numberOfTowerPointsPerTurbine > 0)
            VF_LOG_INFO("tower drag coefficient = {}", this->towerDragCoeff.value());

        if (numberOfHubPointsPerTurbine > 0)
        {
            VF_LOG_INFO("hub drag coefficient = {}", this->hubDragCoeff.value());
            VF_LOG_INFO("hub skin friction coefficient = {}", this->hubSkinFrictionCoeff.value());
        }

    }

    // ActuatorFarmStandalone(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, const real diameter,
    //                        const uint numberOfPointsPerBlade, const std::vector<real>& turbinePositionsX,
    //                        const std::vector<real>& turbinePositionsY, const std::vector<real>& turbinePositionsZ,
    //                        const std::vector<real>& rotorSpeeds, const real smearingWidth, const int level,
    //                        const uint numberOfBlades,
    //                        const std::optional<HubConfig>& hubConfig = std::nullopt,
    //                        const std::optional<TowerConfig>& towerConfig = std::nullopt,
    //                        const std::optional<real> hubDragCoeff = std::nullopt,
    //                        const std::optional<real> hubSkinFrictionCoeff = std::nullopt,
    //                        const std::optional<real> towerDragCoeff = std::nullopt)
    //     : ActuatorFarmStandalone(std::move(para), std::move(cudaMemoryManager), diameter, numberOfPointsPerBlade,
    //                              turbinePositionsX, turbinePositionsY, turbinePositionsZ, rotorSpeeds, smearingWidth,
    //                              level, hubConfig, towerConfig, hubDragCoeff, hubSkinFrictionCoeff, towerDragCoeff,
    //                              numberOfBlades)
    // {
    // }

    ~ActuatorFarmStandalone() override = default;

    void updateForcesAndCoordinates(real time, real deltaT) override;
    static std::vector<real> computeBladeRadii(real diameter, uint numberOfPointsPerBlade);

private:
    void updateBladeForcesAndCoordinates(real deltaT);
    void updateHubForces();
    void updateTowerForces();
    std::vector<real> rotorSpeeds;
    std::optional<real> hubDragCoeff;
    std::optional<real> hubSkinFrictionCoeff;
    std::optional<real> towerDragCoeff;
    std::optional<std::vector<real>> bladeNormalCoefficients;
};

}

#endif

//! \}
