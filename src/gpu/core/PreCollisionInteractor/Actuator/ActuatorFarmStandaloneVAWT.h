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
#ifndef ActuatorFarmStandaloneVAWT_H
#define ActuatorFarmStandaloneVAWT_H

#include "ActuatorFarm.h"
#include <algorithm>
#include <vector>


namespace vf::gpu {

class ActuatorFarmStandaloneVAWT : public ActuatorFarm
{
public:
    ActuatorFarmStandaloneVAWT(
        SPtr<Parameter> para,
        SPtr<CudaMemoryManager> cudaMemoryManager,
        real rotorDiameter,
        uint numberOfBlades,
        uint numberOfPointsPerBlade,
        real rotorHeight,
        const std::vector<real>& positionsTurbineX,
        const std::vector<real>& positionsTurbineY,
        const std::vector<real>& positionsTurbineZ,
        const std::vector<real>& rotorSpeeds,
        real smearingWidth,
        int gridLevelForALM,
        const std::vector<real>& polarAngleOfAttackDeg,
        const std::vector<real>& polarLiftCoefficient,
        const std::vector<real>& polarDragCoefficient,
        real bladeChord,
        real bladePitch,
        real bladeMountingPoint,
        real velocityInlet,
        bool flagLocalizedSmearingWidth = false,
        bool flagFlowCurvature = false,
        bool flagEndEffects = false
    ): ActuatorFarm(std::move(para), std::move(cudaMemoryManager), rotorDiameter,
                 computeBladeRadii(rotorDiameter, numberOfPointsPerBlade), positionsTurbineX, positionsTurbineY,
                 positionsTurbineZ, smearingWidth, gridLevelForALM, true, numberOfBlades, std::nullopt, std::nullopt,
                 ActuatorFarm::VAWTConfig{ rotorHeight, flagLocalizedSmearingWidth }),
    rotorSpeeds(rotorSpeeds),
    rotorHeight(rotorHeight),
    bladeChord(bladeChord),
    bladePitch(bladePitch),
    bladeMountingPoint(bladeMountingPoint),
    velocityInlet(velocityInlet),
    flagFlowCurvature(flagFlowCurvature),
    flagEndEffects(flagEndEffects),
    polarAngleOfAttackDeg(polarAngleOfAttackDeg),
    polarLiftCoefficient(polarLiftCoefficient),
    polarDragCoefficient(polarDragCoefficient),
    bladeHeights(computeBladeHeights(rotorHeight, numberOfPointsPerBlade)),
    endEffectsDistribution(flagEndEffects ? computeEndEffectsDistribution(rotorHeight, bladeChord, numberOfPointsPerBlade)
                        : std::vector<real>(numberOfPointsPerBlade, vf::basics::constant::c1o1))
{
    using namespace vf::basics::constant;

    if (numberOfPointsPerBlade == 0)
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::ActuatorFarmStandaloneVAWT: numberOfPointsPerBlade "
                                 "needs to be > 0.");
    if (numberOfTurbines != this->rotorSpeeds.size())
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::ActuatorFarmStandaloneVAWT: rotor speeds need to have "
                                 "same length as turbine positions!");
    if (this->rotorHeight <= vf::basics::constant::c0o1)
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::ActuatorFarmStandaloneVAWT: rotorHeight needs to be positive!");
    if (this->bladeChord <= vf::basics::constant::c0o1)
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::ActuatorFarmStandaloneVAWT: bladeChord needs to be positive!");
    if (this->polarAngleOfAttackDeg.size() != this->polarLiftCoefficient.size() || this->polarAngleOfAttackDeg.size() != this->polarDragCoefficient.size() || this->polarAngleOfAttackDeg.empty())
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::ActuatorFarmStandaloneVAWT: polarAngleOfAttackDeg/polarLiftCoefficient/polarDragCoefficient vectors need "
                                 "same non-zero size!");
    if (!std::is_sorted(this->polarAngleOfAttackDeg.begin(), this->polarAngleOfAttackDeg.end()))
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::ActuatorFarmStandaloneVAWT: polarAngleOfAttackDeg values need to be sorted "
                                 "in ascending order.");

    VF_LOG_INFO("VAWT parameters:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("rotorHeight [m]    = {}", this->rotorHeight);
    VF_LOG_INFO("bladeChord [m]     = {}", this->bladeChord);
    VF_LOG_INFO("bladePitch [rad]   = {}", this->bladePitch);
    VF_LOG_INFO("mounting point     = {}", this->bladeMountingPoint);
    VF_LOG_INFO("localized smearing = {}", this->flagLocalSmearingWidth ? "on" : "off");
    VF_LOG_INFO("flow curvature     = {}", this->flagFlowCurvature ? "on" : "off");
    VF_LOG_INFO("end effects        = {}", this->flagEndEffects ? "on" : "off");
    this->forceNormal.resize(getTotalNumberOfPoints());
    this->forceTangential.resize(getTotalNumberOfPoints());
    this->angleOfAttackDeg.resize(getTotalNumberOfPoints());
    this->azimuthDeg.resize(getTotalNumberOfPoints());
}

    ~ActuatorFarmStandaloneVAWT() override = default;

    void init() override;
    void updateForcesAndCoordinates(real time, real deltaT) override;

    static std::vector<real> computeBladeRadii(real rotorDiameter, uint numberOfPointsPerBlade);
    static std::vector<real> computeBladeHeights(real rotorHeight, uint numberOfPointsPerBlade);

private:
    static std::vector<real> solveGauss(std::vector<std::vector<real>> D, std::vector<real> b);
    static std::vector<real> computeEndEffectsDistribution(real rotorHeight, real bladeChord, uint numberOfPointsPerBlade);
    static real get_flowCurvature(real vrel, real xchord, real rotorSpeed, real bladeChord);

    void updateCoordinatesVAWT(real time, real deltaT);
    void updateForcesVAWT(real time, real deltaT);
    void appendOutputData(std::vector<std::string>& dataNames, std::vector<std::vector<double>>& nodeData) const override;

    const std::vector<real> rotorSpeeds;
    const real rotorHeight;
    const real bladeChord;
    const real bladePitch;
    const real bladeMountingPoint;
    const real velocityInlet;
    const bool flagFlowCurvature;
    const bool flagEndEffects;
    std::vector<real> forceTangential, forceNormal, angleOfAttackDeg, azimuthDeg;
    const std::vector<real> polarAngleOfAttackDeg;
    const std::vector<real> polarLiftCoefficient;
    const std::vector<real> polarDragCoefficient;
    const std::vector<real> bladeHeights;
    const std::vector<real> endEffectsDistribution;
};

}
#endif

//! \}
