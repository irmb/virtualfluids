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
#ifndef ActuatorFarm_H
#define ActuatorFarm_H

#include "Parameter/Parameter.h"
#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include "logger/Logger.h"
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <cstddef>
#include <optional>
#include <string>
#include <stdexcept>
#include <vector>

namespace vf::gpu {

class GridProvider;

class ActuatorFarm : public PreCollisionInteractor
{
public:
    struct HubConfig {
        real length;
        real radius;
        real positionOffset;
        uint numberOfPointsPerTurbine;
    };
    struct TowerConfig {
        real radius;
        real offset;
        std::vector<real> towerHeights;
        uint numberOfPointsPerTurbine;
    };
    struct VAWTConfig {
        real rotorHeight;
        bool flagLocalSmearingWidth;
    };

    ActuatorFarm(
        SPtr<Parameter> para,
        SPtr<CudaMemoryManager> cudaMemoryManager,
        const real diameter,
        const std::vector<real>& bladeRadii,
        const std::vector<real>& turbinePositionsX,
        const std::vector<real>& turbinePositionsY,
        const std::vector<real>& turbinePositionsZ,
        const real smearingWidth,
        const int level,
        const bool useHostArrays,
        const uint numberOfBlades = 3,
        const std::optional<HubConfig>& hubConfig = std::nullopt,
        const std::optional<TowerConfig>& towerConfig = std::nullopt,
        const std::optional<VAWTConfig>& vawtConfig = std::nullopt
    ) :
        diameter(diameter),
        numberOfBlades(numberOfBlades),
        bladeRadii(bladeRadii),
        numberOfTurbines(static_cast<uint>(turbinePositionsX.size())),
        numberOfPointsPerBlade(static_cast<uint>(bladeRadii.size())),
        numberOfBladePointsPerTurbine(numberOfPointsPerBlade * numberOfBlades),
        initialTurbinePositionsX(turbinePositionsX),
        initialTurbinePositionsY(turbinePositionsY),
        initialTurbinePositionsZ(turbinePositionsZ),
        smearingWidth(smearingWidth),
        level(level),
        useHostArrays(useHostArrays),
        hubLength(hubConfig ? hubConfig->length : 0.0),
        hubRadius(hubConfig ? hubConfig->radius : 0.0),
        hubPositionOffset(hubConfig ? hubConfig->positionOffset : 0.0),
        numberOfHubPointsPerTurbine(hubConfig ? hubConfig->numberOfPointsPerTurbine : 0),
        numberOfHubPoints(this->numberOfHubPointsPerTurbine * numberOfTurbines),
        towerRadius(towerConfig ? towerConfig->radius : 0.0),
        towerOffset(towerConfig ? towerConfig->offset : 0.0),
        numberOfTowerPointsPerTurbine(towerConfig ? towerConfig->numberOfPointsPerTurbine : 0),
        numberOfTowerPoints(this->numberOfTowerPointsPerTurbine * numberOfTurbines),
        towerHeights(towerConfig ? towerConfig->towerHeights : std::vector<real>()),
        numberOfBladePoints(numberOfTurbines * numberOfBladePointsPerTurbine),
        vawtRotorHeight(vawtConfig ? vawtConfig->rotorHeight : 0.0),
        flagLocalSmearingWidth(vawtConfig ? vawtConfig->flagLocalSmearingWidth : false),
        PreCollisionInteractor(std::move(para), std::move(cudaMemoryManager))
    {
        using namespace vf::basics::constant;
        const real deltaX = this->para->getScaledLengthRatio(level);
        const real smearingWidthOverDx = this->smearingWidth / deltaX;
        if(smearingWidthOverDx < c1o1)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: smearing width needs to be larger than dx!");
        if(numberOfTurbines != turbinePositionsY.size() || numberOfTurbines != turbinePositionsZ.size())
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: turbine positions need to have the same length!");
        if(numberOfBlades == 0)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: number of blades must be larger than zero!");
        if(numberOfTowerPointsPerTurbine > 0 && towerHeights.size() < numberOfTurbines)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: tower heights vector size (" + std::to_string(towerHeights.size())
                                     + ") must match number of turbines (" + std::to_string(numberOfTurbines) + ")!");
        if (this->vawtRotorHeight < c0o1)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: VAWT vawtRotorHeight must be non-negative!");

        azimuths = std::vector<real>(numberOfTurbines, c0o1);
        VF_LOG_INFO("ActuatorFarm parameters:");
        VF_LOG_INFO("--------------");
        VF_LOG_INFO("level               = {}", this->level);
        VF_LOG_INFO("number of turbines  = {}", this->numberOfTurbines );
        VF_LOG_INFO("number of blades    = {}", this->numberOfBlades );
        VF_LOG_INFO("points per blade:   = {}", this->numberOfPointsPerBlade);
        VF_LOG_INFO("rotor diameter [m]  = {}", this->diameter);
        VF_LOG_INFO("nodes per diameter  = {}", this->diameter / deltaX);
        VF_LOG_INFO("points per hub:   = {}", this->numberOfHubPointsPerTurbine);
        VF_LOG_INFO("points per tower:   = {}", this->numberOfTowerPointsPerTurbine);
        VF_LOG_INFO("smearing width [m]  = {} ", this->smearingWidth);
        VF_LOG_INFO("smearing width / dx = {} ",this->smearingWidth/deltaX);
        if (this->hasVAWTRotorVolume()) {
            VF_LOG_INFO("VAWT rotor vawtRotorHeight [m] = {}", this->vawtRotorHeight);
            VF_LOG_INFO("VAWT flagLocalSmearingWidth = {}", this->flagLocalSmearingWidth);
        }
    }

    ~ActuatorFarm() override;
    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override;

    void enableOutput(const std::string& outputName, uint tStart, uint tOut) {
        this->outputName = outputName;
        this->writeOutput = true;
        this->tStartOut = tStart;
        this->tOut = tOut;
    }

    void write(const std::string& filename) const;

    uint getNumberOfTurbines() const { return this->numberOfTurbines; };
    uint getNumberOfBladePointsPerTurbine() const { return this->numberOfBladePointsPerTurbine; };
    uint getNumberOfPointsPerBlade() const { return this->numberOfPointsPerBlade; };
    uint getNumberOfBladesPerTurbine() const { return this->numberOfBlades; };

    uint getNumberOfIndices() const { return this->numberOfIndices; };
    uint getNumberOfBladePoints() const { return this->numberOfBladePoints; };
    uint getNumberOfHubPointsPerTurbine() const { return this->numberOfHubPointsPerTurbine; };
    uint getNumberOfTowerPointsPerTurbine() const { return this->numberOfTowerPointsPerTurbine; };
    uint getNumberOfHubPoints() const { return this->numberOfHubPoints; };
    uint getNumberOfTowerPoints() const { return this->numberOfTowerPoints; };

    uint getTotalNumberOfPoints() const {
        return numberOfBladePoints + numberOfHubPoints + numberOfTowerPoints;
    }
    bool hasVAWTRotorVolume() const { return this->vawtRotorHeight > 0.0; }
    bool requiresLocalSmearingWidth() const { return this->hasVAWTRotorVolume() && this->flagLocalSmearingWidth; }
    real getVAWTRotorHeight() const { return this->vawtRotorHeight; }
    bool getFlagLocalSmearingWidth() const { return this->flagLocalSmearingWidth; }

    real* getAllTurbinePosX() const { return turbinePosXH; };
    real* getAllTurbinePosY() const { return turbinePosYH; };
    real* getAllTurbinePosZ() const { return turbinePosZH; };

    real getTurbinePosX(size_t turbine) const { return turbinePosXH[turbine]; };
    real getTurbinePosY(size_t turbine) const { return turbinePosYH[turbine]; };
    real getTurbinePosZ(size_t turbine) const { return turbinePosZH[turbine]; };

    real* getAllBladeCoordsX() const { return this->coordsXH; };
    real* getAllBladeCoordsY() const { return this->coordsYH; };
    real* getAllBladeCoordsZ() const { return this->coordsZH; };
    real* getAllBladeVelocitiesX() const { return this->velocitiesXH; };
    real* getAllBladeVelocitiesY() const { return this->velocitiesYH; };
    real* getAllBladeVelocitiesZ() const { return this->velocitiesZH; };
    real* getAllBladeForcesX() const { return this->forcesXH; };
    real* getAllBladeForcesY() const { return this->forcesYH; };
    real* getAllBladeForcesZ() const { return this->forcesZH; };
    real* getAllBladeSmearingWidth() const { return this->localSmearingWidthH; };
    real* getAllBladeLocalSmearingWidth() const { return this->getAllBladeSmearingWidth(); };

    real* getAllHubCoordsX() const { return numberOfHubPoints > 0 ? &this->coordsXH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubCoordsY() const { return numberOfHubPoints > 0 ? &this->coordsYH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubCoordsZ() const { return numberOfHubPoints > 0 ? &this->coordsZH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubVelocitiesX() const { return numberOfHubPoints > 0 ? &this->velocitiesXH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubVelocitiesY() const { return numberOfHubPoints > 0 ? &this->velocitiesYH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubVelocitiesZ() const { return numberOfHubPoints > 0 ? &this->velocitiesZH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubForcesX() const { return numberOfHubPoints > 0 ? &this->forcesXH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubForcesY() const { return numberOfHubPoints > 0 ? &this->forcesYH[this->numberOfBladePoints] : nullptr; };
    real* getAllHubForcesZ() const { return numberOfHubPoints > 0 ? &this->forcesZH[this->numberOfBladePoints] : nullptr; };

    real* getAllTowerCoordsX() const { return numberOfTowerPoints > 0 ? &this->coordsXH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerCoordsY() const { return numberOfTowerPoints > 0 ? &this->coordsYH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerCoordsZ() const { return numberOfTowerPoints > 0 ? &this->coordsZH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerVelocitiesX() const { return numberOfTowerPoints > 0 ? &this->velocitiesXH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerVelocitiesY() const { return numberOfTowerPoints > 0 ? &this->velocitiesYH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerVelocitiesZ() const { return numberOfTowerPoints > 0 ? &this->velocitiesZH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerForcesX() const { return numberOfTowerPoints > 0 ? &this->forcesXH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerForcesY() const { return numberOfTowerPoints > 0 ? &this->forcesYH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };
    real* getAllTowerForcesZ() const { return numberOfTowerPoints > 0 ? &this->forcesZH[this->numberOfBladePoints+this->numberOfHubPoints] : nullptr; };

    real* getAllBladeCoordsXDevice() const { return this->coordsXDCurrentTimestep; };
    real* getAllBladeCoordsYDevice() const { return this->coordsYDCurrentTimestep; };
    real* getAllBladeCoordsZDevice() const { return this->coordsZDCurrentTimestep; };
    real* getAllBladeVelocitiesXDevice() const { return this->velocitiesXDCurrentTimestep; };
    real* getAllBladeVelocitiesYDevice() const { return this->velocitiesYDCurrentTimestep; };
    real* getAllBladeVelocitiesZDevice() const { return this->velocitiesZDCurrentTimestep; };
    real* getAllBladeForcesXDevice() const { return this->forcesXDCurrentTimestep; };
    real* getAllBladeForcesYDevice() const { return this->forcesYDCurrentTimestep; };
    real* getAllBladeForcesZDevice() const { return this->forcesZDCurrentTimestep; };
    real* getAllBladeSmearingWidthDevice() const { return this->localSmearingWidthDCurrentTimestep; };
    real* getAllBladeLocalSmearingWidthDevice() const { return this->getAllBladeSmearingWidthDevice(); };

    real* getAllHubCoordsXDevice() const { return numberOfHubPoints > 0 ? &this->coordsXDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubCoordsYDevice() const { return numberOfHubPoints > 0 ? &this->coordsYDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubCoordsZDevice() const { return numberOfHubPoints > 0 ? &this->coordsZDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubVelocitiesXDevice() const { return numberOfHubPoints > 0 ? &this->velocitiesXDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubVelocitiesYDevice() const { return numberOfHubPoints > 0 ? &this->velocitiesYDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubVelocitiesZDevice() const { return numberOfHubPoints > 0 ? &this->velocitiesZDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubForcesXDevice() const { return numberOfHubPoints > 0 ? &this->forcesXDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubForcesYDevice() const { return numberOfHubPoints > 0 ? &this->forcesYDCurrentTimestep[numberOfBladePoints] : nullptr; };
    real* getAllHubForcesZDevice() const { return numberOfHubPoints > 0 ? &this->forcesZDCurrentTimestep[numberOfBladePoints] : nullptr; };

    real* getAllTowerCoordsXDevice() const { return numberOfTowerPoints > 0 ? &this->coordsXDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerCoordsYDevice() const { return numberOfTowerPoints > 0 ? &this->coordsYDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerCoordsZDevice() const { return numberOfTowerPoints > 0 ? &this->coordsZDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerVelocitiesXDevice() const { return numberOfTowerPoints > 0 ? &this->velocitiesXDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerVelocitiesYDevice() const { return numberOfTowerPoints > 0 ? &this->velocitiesYDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerVelocitiesZDevice() const { return numberOfTowerPoints > 0 ? &this->velocitiesZDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerForcesXDevice() const { return numberOfTowerPoints > 0 ? &this->forcesXDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerForcesYDevice() const { return numberOfTowerPoints > 0 ? &this->forcesYDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };
    real* getAllTowerForcesZDevice() const { return numberOfTowerPoints > 0 ? &this->forcesZDCurrentTimestep[numberOfBladePoints+numberOfHubPoints] : nullptr; };

    void setAllBladeCoords(const real* bladeCoordsX, const real* bladeCoordsY, const real* bladeCoordsZ) const;
    void setAllBladeVelocities(const real* bladeVelocitiesX, const real* bladeVelocitiesY, const real* bladeVelocitiesZ) const;
    void setAllBladeForces(const real* bladeForcesX, const real* bladeForcesY, const real* bladeForcesZ) const;

    void setTurbineBladeCoords(size_t turbine, const real* bladeCoordsX, const real* bladeCoordsY, const real* bladeCoordsZ) const;
    void setTurbineBladeVelocities(size_t turbine, const real* bladeVelocitiesX, const real* bladeVelocitiesY, const real* bladeVelocitiesZ) const;
    void setTurbineBladeForces(size_t turbine, const real* bladeForcesX, const real* bladeForcesY, const real* bladeForcesZ) const;

    void setTurbineAzimuth(size_t turbine, real azimuth){azimuths[turbine] = azimuth;}

    virtual void updateForcesAndCoordinates(real time, real deltaT)=0;
    virtual void appendOutputData(std::vector<std::string>&,
                                  std::vector<std::vector<double>>&) const {};

private:
    void initTurbineGeometries();
    void initBoundingVolumes();
    void initCoords();
    void initVelocities();
    void initForces();
    void initIndices();
    std::string getFilename(uint t) const;
    void swapDeviceArrays();
    void generateHubAxisPoints(uint turbineIndex, uint& pointIndex);
    void generateTowerAxisPoints(uint turbineIndex, uint& pointIndex);
public:
    real* coordsXH, *coordsYH, *coordsZH;
    real* velocitiesXH, *velocitiesYH, *velocitiesZH;
    real* forcesXH, *forcesYH, *forcesZH;
    real* localSmearingWidthH{nullptr};
    uint* indicesH;

    uint* boundingVolumeIndicesH;
    real* turbinePosXH, *turbinePosYH, *turbinePosZH;

    real* coordsXDCurrentTimestep, * coordsYDCurrentTimestep, * coordsZDCurrentTimestep;
    real* coordsXDPreviousTimestep, * coordsYDPreviousTimestep, * coordsZDPreviousTimestep;
    real* velocitiesXDCurrentTimestep, * velocitiesYDCurrentTimestep, * velocitiesZDCurrentTimestep;
    real* velocitiesXDPreviousTimestep, * velocitiesYDPreviousTimestep, * velocitiesZDPreviousTimestep;
    real* forcesXDCurrentTimestep, * forcesYDCurrentTimestep, * forcesZDCurrentTimestep;
    real* forcesXDPreviousTimestep, * forcesYDPreviousTimestep, * forcesZDPreviousTimestep;
    real* localSmearingWidthDCurrentTimestep{nullptr}, *localSmearingWidthDPreviousTimestep{nullptr};
    uint* indicesD;

    real* turbinePosXD, *turbinePosYD, *turbinePosZD;
    uint* boundingVolumeIndicesD;

protected:
    real getRotorBoundingSmearingWidth() const
    {
        return (this->useLocalSmearingWidth() && this->vawtBoundingSmearingWidth > 0.0)
            ? this->vawtBoundingSmearingWidth
            : this->smearingWidth;
    }
    real getVAWTRotorBoundingMargin() const
    {
        using namespace vf::basics::constant;
        return this->useLocalSmearingWidth()
            ? c3o1 * this->getRotorBoundingSmearingWidth()
            : c3o1 * this->smearingWidth;
    }
    bool useLocalSmearingWidth() const
    {
        return this->hasVAWTRotorVolume() && this->flagLocalSmearingWidth;
    }

    std::vector<real> bladeRadii, initialTurbinePositionsX, initialTurbinePositionsY, initialTurbinePositionsZ;
    std::vector<real> azimuths;
    const real diameter;
    const bool useHostArrays;
    const uint numberOfTurbines, numberOfPointsPerBlade, numberOfBlades, numberOfBladePointsPerTurbine;
    const real smearingWidth; // in m
    const int level;
    uint numberOfIndices{0};
    const uint numberOfBladePoints;
    const real vawtRotorHeight;
    const bool flagLocalSmearingWidth;
    real vawtBoundingSmearingWidth{0.0};
    const real hubLength, hubRadius, hubPositionOffset;
    const uint numberOfHubPointsPerTurbine;
    const uint numberOfHubPoints;
    const real towerRadius, towerOffset;
    const uint numberOfTowerPointsPerTurbine;
    const uint numberOfTowerPoints;
    const std::vector<real> towerHeights;
    real maxTowerHeight{0.0};
    int streamIndex;

    bool writeOutput{false};
    std::string outputName;
    uint tOut{0};
    uint tStartOut{0};
};

}

#endif

//! \}
