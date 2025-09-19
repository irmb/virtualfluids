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

#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <stdexcept>
#include <vector>

class GridProvider;

class ActuatorFarm : public PreCollisionInteractor
{
public:
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
        const real deltaT,
        const real deltaX,
        const bool useHostArrays
    ) :
        diameter(diameter),
        bladeRadii(bladeRadii),
        numberOfNodesPerBlade(static_cast<uint>(bladeRadii.size())),
        numberOfNodesPerTurbine(numberOfNodesPerBlade*numberOfBlades),
        numberOfTurbines(static_cast<uint>(turbinePositionsX.size())),
        initialTurbinePositionsX(turbinePositionsX),
        initialTurbinePositionsY(turbinePositionsY),
        initialTurbinePositionsZ(turbinePositionsZ),
        smearingWidth(smearingWidth),
        level(level),
        useHostArrays(useHostArrays),
        deltaT(deltaT*std::exp2(-level)),
        deltaX(deltaX*std::exp2(-level)),
        invDeltaX(vf::basics::constant::c1o1/deltaX),
        PreCollisionInteractor(std::move(para), std::move(cudaMemoryManager))
    {
        if(this->smearingWidth < this->deltaX)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: smearing width needs to be larger than dx!");
        if(numberOfTurbines != turbinePositionsY.size() || numberOfTurbines != turbinePositionsZ.size())
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: turbine positions need to have the same length!");
        azimuths = std::vector<real>(numberOfTurbines, 0.0);
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

    real getDeltaT() const { return this->deltaT; };
    real getDeltaX() const { return this->deltaX; };

    uint getNumberOfTurbines() const { return this->numberOfTurbines; };
    uint getNumberOfNodesPerTurbine() const { return this->numberOfNodesPerTurbine; };
    uint getNumberOfNodesPerBlade() const { return this->numberOfNodesPerBlade; };
    uint getNumberOfBladesPerTurbine() const { return ActuatorFarm::numberOfBlades; };

    uint getNumberOfIndices() const { return this->numberOfIndices; };
    uint getNumberOfGridNodes() const { return this->numberOfGridNodes; };

    real* getAllTurbinePosX() const { return turbinePosXH; };
    real* getAllTurbinePosY() const { return turbinePosYH; };
    real* getAllTurbinePosZ() const { return turbinePosZH; };

    real getTurbinePosX(size_t turbine) const { return turbinePosXH[turbine]; };
    real getTurbinePosY(size_t turbine) const { return turbinePosYH[turbine]; };
    real getTurbinePosZ(size_t turbine) const { return turbinePosZH[turbine]; };

    real* getAllBladeCoordsX() const { return this->bladeCoordsXH; };
    real* getAllBladeCoordsY() const { return this->bladeCoordsYH; };
    real* getAllBladeCoordsZ() const { return this->bladeCoordsZH; };
    real* getAllBladeVelocitiesX() const { return this->bladeVelocitiesXH; };
    real* getAllBladeVelocitiesY() const { return this->bladeVelocitiesYH; };
    real* getAllBladeVelocitiesZ() const { return this->bladeVelocitiesZH; };
    real* getAllBladeForcesX() const { return this->bladeForcesXH; };
    real* getAllBladeForcesY() const { return this->bladeForcesYH; };
    real* getAllBladeForcesZ() const { return this->bladeForcesZH; };

    real* getTurbineBladeCoordsX(size_t turbine) const { return &this->bladeCoordsXH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsY(size_t turbine) const { return &this->bladeCoordsYH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsZ(size_t turbine) const { return &this->bladeCoordsZH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesX(size_t turbine) const { return &this->bladeVelocitiesXH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesY(size_t turbine) const { return &this->bladeVelocitiesYH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesZ(size_t turbine) const { return &this->bladeVelocitiesZH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesX(size_t turbine) const { return &this->bladeForcesXH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesY(size_t turbine) const { return &this->bladeForcesYH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesZ(size_t turbine) const { return &this->bladeForcesZH[turbine*numberOfNodesPerTurbine]; };

    real* getAllBladeCoordsXDevice() const { return this->bladeCoordsXDCurrentTimestep; };
    real* getAllBladeCoordsYDevice() const { return this->bladeCoordsYDCurrentTimestep; };
    real* getAllBladeCoordsZDevice() const { return this->bladeCoordsZDCurrentTimestep; };
    real* getAllBladeVelocitiesXDevice() const { return this->bladeVelocitiesXDCurrentTimestep; };
    real* getAllBladeVelocitiesYDevice() const { return this->bladeVelocitiesYDCurrentTimestep; };
    real* getAllBladeVelocitiesZDevice() const { return this->bladeVelocitiesZDCurrentTimestep; };
    real* getAllBladeForcesXDevice() const { return this->bladeForcesXDCurrentTimestep; };
    real* getAllBladeForcesYDevice() const { return this->bladeForcesYDCurrentTimestep; };
    real* getAllBladeForcesZDevice() const { return this->bladeForcesZDCurrentTimestep; };

    real* getTurbineBladeCoordsXDevice(size_t turbine) const { return &this->bladeCoordsXDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsYDevice(size_t turbine) const { return &this->bladeCoordsYDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsZDevice(size_t turbine) const { return &this->bladeCoordsZDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesXDevice(size_t turbine) const { return &this->bladeVelocitiesXDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesYDevice(size_t turbine) const { return &this->bladeVelocitiesYDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesZDevice(size_t turbine) const { return &this->bladeVelocitiesZDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesXDevice(size_t turbine) const { return &this->bladeForcesXDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesYDevice(size_t turbine) const { return &this->bladeForcesYDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesZDevice(size_t turbine) const { return &this->bladeForcesZDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };

    void setAllBladeCoords(const real* bladeCoordsX, const real* bladeCoordsY, const real* bladeCoordsZ) const;
    void setAllBladeVelocities(const real* bladeVelocitiesX, const real* bladeVelocitiesY, const real* bladeVelocitiesZ) const;
    void setAllBladeForces(const real* bladeForcesX, const real* bladeForcesY, const real* bladeForcesZ) const;

    void setTurbineBladeCoords(size_t turbine, const real* bladeCoordsX, const real* bladeCoordsY, const real* bladeCoordsZ) const;
    void setTurbineBladeVelocities(size_t turbine, const real* bladeVelocitiesX, const real* bladeVelocitiesY, const real* bladeVelocitiesZ) const;
    void setTurbineBladeForces(size_t turbine, const real* bladeForcesX, const real* bladeForcesY, const real* bladeForcesZ) const;

    void setTurbineAzimuth(size_t turbine, real azimuth){azimuths[turbine] = azimuth;}

    virtual void updateForcesAndCoordinates()=0;

private:
    void initTurbineGeometries();
    void initBoundingSpheres();
    void initBladeCoords();
    void initBladeVelocities();
    void initBladeForces();
    void initBladeIndices();
    std::string getFilename(uint t) const;
    void swapDeviceArrays();

public:
    real* bladeCoordsXH, * bladeCoordsYH, * bladeCoordsZH;
    real* bladeCoordsXDPreviousTimestep, * bladeCoordsYDPreviousTimestep, * bladeCoordsZDPreviousTimestep;
    real* bladeCoordsXDCurrentTimestep, * bladeCoordsYDCurrentTimestep, * bladeCoordsZDCurrentTimestep;    
    real* bladeVelocitiesXH, * bladeVelocitiesYH, * bladeVelocitiesZH;
    real* bladeVelocitiesXDPreviousTimestep, * bladeVelocitiesYDPreviousTimestep, * bladeVelocitiesZDPreviousTimestep;
    real* bladeVelocitiesXDCurrentTimestep, * bladeVelocitiesYDCurrentTimestep, * bladeVelocitiesZDCurrentTimestep;
    real* bladeForcesXH, * bladeForcesYH, * bladeForcesZH;
    real* bladeForcesXDPreviousTimestep, * bladeForcesYDPreviousTimestep, * bladeForcesZDPreviousTimestep;
    real* bladeForcesXDCurrentTimestep, * bladeForcesYDCurrentTimestep, * bladeForcesZDCurrentTimestep;
    uint* bladeIndicesH;
    uint* bladeIndicesD; 
    uint* boundingSphereIndicesH;
    uint* boundingSphereIndicesD;
    real* turbinePosXH, *turbinePosYH, *turbinePosZH;
    real* turbinePosXD, *turbinePosYD, *turbinePosZD;

protected:
    static constexpr uint numberOfBlades{3};
    std::vector<real> bladeRadii, initialTurbinePositionsX, initialTurbinePositionsY, initialTurbinePositionsZ;
    std::vector<real> azimuths;
    const real diameter;
    const bool useHostArrays;
    const real deltaT, deltaX, invDeltaX;
    const uint numberOfTurbines, numberOfNodesPerBlade, numberOfNodesPerTurbine;
    const real smearingWidth; // in m
    const int level;
    uint numberOfIndices{0};
    uint numberOfGridNodes{0};

    real factorGaussian;
    int streamIndex;

    bool writeOutput{false};
    std::string outputName;
    uint tOut{0};
    uint tStartOut{0};
};

#endif

//! \}
