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
//! \author Henrik Asmuth, Henry Korb
//=======================================================================================
#define _USE_MATH_DEFINES
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <new>
#include <random>
#include <stdexcept>
#include <string>

//////////////////////////////////////////////////////////////////////////

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>
#include <basics/geometry3d/Axis.h>

#include <logger/Logger.h>

#include <gpu/core/BoundaryConditions/BoundaryConditionFactory.h>
#include <gpu/core/Calculation/Simulation.h>
#include <gpu/core/Cuda/CudaMemoryManager.h>
#include <gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h>
#include <gpu/core/GridScaling/GridScalingFactory.h>
#include <gpu/core/Kernel/KernelTypes.h>
#include <gpu/core/Parameter/Parameter.h>
#include <gpu/core/PreCollisionInteractor/BuoyancyProvider/BuoyancyProvider.h>
#include <gpu/core/PreCollisionInteractor/CoriolisForce.h>
#include <gpu/core/PreCollisionInteractor/DampingLayer/DampingLayer.h>
#include <gpu/core/Samplers/PlanarAverageProbe.h>
#include <gpu/core/Samplers/Probe.h>
#include <gpu/core/Samplers/WallModelProbe.h>
#include <gpu/core/TurbulenceModels/TurbulenceModelFactory.h>

#include <gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h>
#include <gpu/GridGenerator/grid/BoundaryConditions/Side.h>
#include <gpu/GridGenerator/grid/GridBuilder/MultipleGridBuilder.h>
#include <gpu/GridGenerator/grid/MultipleGridBuilderFacade.h>

#include <parallel/MPICommunicator.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace vf::basics::constant;
using temperatureProfileFunc = std::function<real(real)>;
using velocityProfileFunc = std::function<std::tuple<real, real, real, real>(real, real, real)>;

const std::string defaultConfigFile("sabl.cfg");
const std::string simulationName("BoundaryLayer");

temperatureProfileFunc getTemperatureFuncConstantLapseRate(real surfaceTemperature, real lapseRate)
{
    return [=](real coordZ) { return surfaceTemperature + coordZ * lapseRate; };
}

temperatureProfileFunc getTemperatureFuncTwoLayers(real surfaceTemperature, real inversionHeight, real lapseRate)
{
    return [=](real coordZ) { return surfaceTemperature + lapseRate * std::max(coordZ - inversionHeight, c0o1); };
}

//! follows <a href="https://doi.org/10.1175/1520-0450(2004)043<0925:AMTDTC>2.0.CO;2"> Rampanelli and Zardi (2004) </a>
temperatureProfileFunc getTemperatureFuncAllaerts(real inversionHeight, real inversionThickness, real inversionStrength,
                                                  real lapseRate, real surfaceTemperature)
{
    return [=](real coordZ) {
        const real eta = c3o1 * (coordZ - inversionHeight) / inversionThickness;
        const real fBelowInversion = (std::tanh(eta) + c1o1) * c1o2;
        const real fAboveInversion = (std::log(c2o1 * std::cosh(eta)) + eta) * c1o2;
        return inversionStrength * fBelowInversion + c1o3 * lapseRate * inversionThickness * fAboveInversion +
               surfaceTemperature;
    };
}

velocityProfileFunc getConstantVelocityProfile(real geostrophicWindX, real geostrophicWindY)
{
    return [=](real, real, real) { return std::tuple { geostrophicWindX, geostrophicWindY, c0o1, c0o1 }; };
}

velocityProfileFunc getVelocityProfileFuncAllaerts(real frictionVelocity, real roughnessLength, real rotationAngle,
                                                   real inversionHeight, real coriolisParameter, real mergingRegionThickness,
                                                   real geostrophicWind, real vonKarmanConstant)
{
    return [=](real, real, real coordZ) {
        const real zeta = coordZ / inversionHeight;

        const real f_u = 1.57F * zeta - 2.68F * zeta * zeta;
        const real f_v = 13.2F * zeta - 8.7F * zeta * zeta;
        const real u_ABL = frictionVelocity / vonKarmanConstant * (std::log1p(coordZ / roughnessLength) + f_u);
        const real v_ABL = frictionVelocity / vonKarmanConstant * f_v * (coriolisParameter >= c0o1 ? c1o1 : -c1o1);

        const real merging_factor = std::tanh((zeta - c1o2) * c2o1 * inversionHeight / mergingRegionThickness);

        const real u = c1o2 * (u_ABL * (c1o1 - merging_factor) + geostrophicWind * (c1o1 + merging_factor));
        const real v = c1o2 * (v_ABL * (c1o1 - merging_factor) + c0o1);

        return std::tuple { u * std::cos(rotationAngle) - v * std::sin(rotationAngle),
                            u * std::sin(rotationAngle) + v * std::cos(rotationAngle), c0o1, c0o1 };
    };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void run(const vf::basics::ConfigurationFile& config)
{
    const real gravity = 9.81F;
    const real vonKarmanConstant = 0.4F;
    const bool useCoriolis = true;
    const bool useSurfaceLayer = true;
    const real mergingRegionThickness = c2o1 * c100o1;

    const real lengthX = config.getValue("L_x", 6000.0);
    const real lengthY = config.getValue("L_y", 3000.0);
    const real lengthZ = config.getValue("L_z", 1000.0);
    const real boundaryLayerHeight = config.getValue("BoundaryLayerHeight", 1000.0);

    const real viscosity = config.getValue("viscosity", 1.56e-5F);
    const real mach = config.getValue<real>("Ma", 0.1);
    const uint nodesPerBoundaryLayerHeight = config.getValue<uint>("nz", 64);

    const real frictionVelocity = config.getValue("FrictionVelocity", 0.28F);
    const real roughnessLength = config.getValue("RoughnessLength", 2.0e-4F);
    const uint samplingOffset = config.getValue<uint>("SamplingOffset", 3);
    const real prandtlNumber = config.getValue<real>("PrandtlNumber", 0.7F);
    const real referenceTemperature = config.getValue<real>("ReferenceTemperature", 288);
    const real surfaceTemperature = config.getValue<real>("SurfaceTemperature", 288) - referenceTemperature;
    const real inversionStrength = config.getValue<real>("InversionStrength", 2.5);
    const real inversionHeight = config.getValue<real>("InversionHeight", c1o2 * boundaryLayerHeight);
    const real inversionThickness = config.getValue<real>("InversionThickness", 100.0);
    const real coriolisParameter = config.getValue<real>("CoriolisParameter", 1e-4F);
    const real lapseRate = config.getValue<real>("LapseRate", 1.0) * 1e-3F;               // is conventionally given in K/km
    const real heatingRate = config.getValue<real>("HeatingRate", c0o1) * c1o36 * c1o100; // is conventionally given in K/h
    const real surfaceHeatFlux = config.getValue<real>("SurfaceHeatFlux", c0o1);
    const bool useHeatFlux = config.getValue<bool>("UseHeatFlux", false);
    const std::string initialVelocityProfile = config.getValue<std::string>("InitialVelocityProfile", "Constant");
    const std::string initialTemperatureProfile =
        config.getValue<std::string>("InitialTemperatureProfile", "ConstantLapseRate");
    const real geostrophicWind = config.getValue<real>("GeostrophicWindSpeed", 10.0);
    const real windDirection = config.getValue<real>("GeostrophicWindDirection", 0.0) * cPio180;
    const real temperaturePerturbationHeight = config.getValue<real>("TemperaturePerturbationHeight", 200);
    const real temperaturePerturbationStrength = config.getValue<real>("TemperaturePerturbationStrength", 0.1);

    // all in s
    const real tStartOut = config.getValue<real>("tStartOut");
    const real tOut = config.getValue<real>("tOut");
    const real tEnd = config.getValue<real>("tEnd"); // total time of simulation

    const real tStartSampling = config.getValue<real>("tStartSampling");
    const real tStartTmpAveraging = config.getValue<real>("tStartTmpAveraging");
    const real tSampling = config.getValue<real>("tSampling");
    const real tStartOutProbe = config.getValue<real>("tStartOutProbe");
    const real tOutProbe = config.getValue<real>("tOutProbe");

    const real dampingLayerStartHeight = config.getValue<real>("DampingLayerStartHeight", 0.7 * lengthZ);
    const real dampingFactor = config.getValue<real>("DampingFactor", 1e6);
    const bool useRayleighDamping = config.getValue("useRayleighDamping", false);

    const real velocity = geostrophicWind;

    const real geostrophicWindX = geostrophicWind * std::cos(windDirection);
    const real geostrophicWindY = geostrophicWind * std::sin(windDirection);

    const real deltaX = boundaryLayerHeight / real(nodesPerBoundaryLayerHeight);
    const real deltaT = deltaX * mach / (std::sqrt(3) * velocity);
    const real velocityLB = velocity * deltaT / deltaX;              // LB units
    const real viscosityLB = viscosity * deltaT / (deltaX * deltaX); // LB units
    const real diffusivity = viscosity / prandtlNumber;
    const real diffusivityLB = diffusivity * deltaT / (deltaX * deltaX);
    const real bruntVaisalaFrequency = std::sqrt(gravity / referenceTemperature * lapseRate);

    temperatureProfileFunc initialTemperatureProfileFunc;
    velocityProfileFunc initialVelocityProfileFunc;
    if (initialTemperatureProfile == "ConstantLapseRate") {
        initialTemperatureProfileFunc = getTemperatureFuncConstantLapseRate(surfaceTemperature, lapseRate);
    } else if (initialTemperatureProfile == "TwoLayers") {
        initialTemperatureProfileFunc = getTemperatureFuncTwoLayers(surfaceTemperature, inversionHeight, lapseRate);
    } else if (initialTemperatureProfile == "ALLAERTS") {
        initialTemperatureProfileFunc = getTemperatureFuncAllaerts(inversionHeight, inversionThickness, inversionStrength,
                                                                   lapseRate, surfaceTemperature);
    } else
        throw std::runtime_error(
            "Unrecognized temperature initialization profile! Must be TwoLayers, ALLAERTS or ConstantLapseRate");

    if (initialVelocityProfile == "Constant")
        initialVelocityProfileFunc = getConstantVelocityProfile(geostrophicWindX, geostrophicWindY);
    else if (initialVelocityProfile == "ALLAERTS")
        initialVelocityProfileFunc =
            getVelocityProfileFuncAllaerts(frictionVelocity, roughnessLength, windDirection, inversionHeight,
                                           coriolisParameter, mergingRegionThickness, geostrophicWind, vonKarmanConstant);
    else
        throw std::runtime_error("Unrecognized velocity initialization profile, must be Constant or ALLAERTS");

    const real temperatureTop = initialTemperatureProfileFunc(lengthZ);

    VF_LOG_INFO("velocity  [dx/dt] = {:3.3f}", velocityLB);
    VF_LOG_INFO("dt   = {:3.2f}s", deltaT);
    VF_LOG_INFO("dx   = {:3.1f}m", deltaX);
    VF_LOG_INFO("viscosity [dx^2/dt] = {:3.2e}", viscosityLB);
    VF_LOG_INFO("diffusivityLB = {:3.2e}", diffusivityLB);
    VF_LOG_INFO("u* /(dx/dt) = {:3.2e}", frictionVelocity * deltaT / deltaX);
    VF_LOG_INFO("Reference Temperature = {:3.2f}K", referenceTemperature);
    VF_LOG_INFO("Surface Temperature = {:3.2f}K", surfaceTemperature + referenceTemperature);
    VF_LOG_INFO("z0 = {:3.3e}m", roughnessLength);
    VF_LOG_INFO("heating rate = {:3.3e}K/s", heatingRate);
    VF_LOG_INFO("inversionHeight = {:4.0f}m", inversionHeight);
    VF_LOG_INFO("inversionThickness = {:3.1f}m", inversionThickness);
    VF_LOG_INFO("inversionStrength = {:3.2f}K", inversionStrength);
    VF_LOG_INFO("geostrophicWind = {:3.2f}, {:3.2f} m/s", geostrophicWindX, geostrophicWindY);
    if (useCoriolis)
        VF_LOG_INFO("coriolisParameter = {:3.2e}s^-1", coriolisParameter);
    VF_LOG_INFO("lapseRate = {:3.2f}K/km", lapseRate * 1e3);
    VF_LOG_INFO("Brunt-Väisälä Frequency = {:3.2e}s^-1", bruntVaisalaFrequency);
    VF_LOG_INFO("T_BC = {:3.2f}", temperatureTop);
    if (useHeatFlux)
        VF_LOG_INFO("Use constant heat flux of {:3.2e}", surfaceHeatFlux);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Grid and boundary conditions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto& communicator = *vf::parallel::MPICommunicator::getInstance();
    const uint numberOfProcesses = communicator.getNumberOfProcesses();
    const uint processID = communicator.getProcessID();

    auto gridDimensions = std::make_shared<GridDimensions>(c0o1, lengthX, c0o1, lengthY, c0o1, lengthZ, deltaX);
    auto gridBuilderFacade = std::make_unique<MultipleGridBuilderFacade>(gridDimensions, c1o1 * deltaX);

    const real subdomainLength = lengthX / static_cast<real>(numberOfProcesses);
    for (uint subdomain = 1; subdomain < numberOfProcesses; subdomain++)
        gridBuilderFacade->addDomainSplit(subdomain * subdomainLength, Axis::x);

    gridBuilderFacade->setPeriodicBoundaryCondition(true, true, false);
    // gridBuilderFacade->setPeriodicShiftOnXBoundaryInYDirection((c1o10 + c1o1) * lengthZ);
    gridBuilderFacade->createGrids(processID);

    if (!useSurfaceLayer) {
        gridBuilderFacade->setADFluxBoundaryCondition(SideType::MZ, c0o1, c0o1, c0o1, c0o1, deltaX);
        gridBuilderFacade->setStressBoundaryCondition(SideType::MZ, c0o1, c0o1, c1o1, samplingOffset, vonKarmanConstant,
                                                      roughnessLength, deltaX);
    } else {
        gridBuilderFacade->setSurfaceLayerBoundaryCondition(SideType::MZ, c0o1, c0o1, c1o1, samplingOffset, vonKarmanConstant,
                                                            roughnessLength, roughnessLength, surfaceHeatFlux,
                                                            surfaceTemperature, heatingRate, deltaX, deltaT);
    }

    gridBuilderFacade->setSlipBoundaryCondition(SideType::PZ, c0o1, c0o1, -c1o1);
    gridBuilderFacade->setADFluxBoundaryCondition(SideType::PZ, c0o1, c0o1, -c1o1, -lapseRate, deltaX);

    BoundaryConditionFactory bcFactory;
    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackCompressible);
    const auto fluxBC = useHeatFlux ? BoundaryConditionFactory::SurfaceLayerBC::SurfaceHeatFlux
                                    : BoundaryConditionFactory::SurfaceLayerBC::SurfaceTemperature;
    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackCompressible, fluxBC);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipBounceBack);
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityBounceBack);
    bcFactory.setAdvectionDiffusionFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxBounceBack);

    auto para = std::make_shared<Parameter>(numberOfProcesses, processID, &config);
    para->worldLength = lengthX + deltaX;
    para->setMaxDev(10);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initial Conditions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz) {
        auto [vxSI, vySI, vzSI, rhoSI] = initialVelocityProfileFunc(coordX, coordY, coordZ);
        vx = vxSI * deltaT / deltaX;
        vy = vySI * deltaT / deltaX;
        vz = vzSI * deltaT / deltaX;
        rho = rhoSI;
    });

    auto randomDistributionTemp =
        std::uniform_real_distribution<real>(-temperaturePerturbationStrength, temperaturePerturbationStrength);
    auto randomGenerator = std::default_random_engine(42);

    para->setInitialConditionAD([&](real, real, real coordZ) {
        if (coordZ < temperaturePerturbationHeight)
            return initialTemperatureProfileFunc(coordZ) + randomDistributionTemp(randomGenerator);
        return initialTemperatureProfileFunc(coordZ);
    });

    para->setInitialLocalReferenceTemperature(
        [&](real, real, real coordZ) -> real { return initialTemperatureProfileFunc(coordZ); });

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set parameters
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);
    para->setUseStreams(true);

    para->setForcing(0, 0, 0);
    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(deltaX / deltaT);
    para->setViscosityRatio(deltaX * deltaX / deltaT);
    para->setDensityRatio(1.0);

    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);

    para->setTimestepStartOut(uint(tStartOut / deltaT));
    para->setTimestepOut(uint(tOut / deltaT));
    para->setTimestepEnd(uint(tEnd / deltaT));
    para->setGravity(gravity * deltaT * deltaT / deltaX);
    para->setBuoyancyFactor(para->getGravity() / referenceTemperature);

    // Advection Diffusion
    para->setDiffOn(true);
    para->setADKernel(vf::advectionDiffusionKernel::compressible::F16);
    para->setDiffusivity(diffusivityLB);
    para->setIsBodyForce(true);
    para->setBuoyancyEnabled(true);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Add Actuators
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    para->addInteractor(std::make_shared<BuoyancyProviderPlanarAverage>(para, cudaMemoryManager));

    if (useCoriolis)
        para->addInteractor(
            std::make_shared<CoriolisForce>(para, cudaMemoryManager, geostrophicWindX, geostrophicWindY, coriolisParameter));

    if (useRayleighDamping) {
        const auto dampingFunction = [&](real zNormalized) -> real {
            return dampingFactor * deltaT * std::sin(cPi * c1o2 * zNormalized) * std::sin(cPi * c1o2 * zNormalized);
        };
        para->addInteractor(std::make_shared<DampingLayer>(para, cudaMemoryManager, DampingLayer::DampingLayerType::Rayleigh,
                                                           Axis::z, dampingFunction, dampingLayerStartHeight,
                                                           lengthZ + deltaX));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Add Probes
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto planarAverageProbe = std::make_shared<PlanarAverageProbe>(
        para, cudaMemoryManager, para->getOutputPath(), "planeProbe", tStartSampling / deltaT, tStartTmpAveraging / deltaT,
        tSampling / deltaT, tStartOutProbe / deltaT, tOutProbe / deltaT, Axis::z, true, true);
    planarAverageProbe->addAllAvailableStatistics();
    planarAverageProbe->setFileNameToNOut();
    para->addSampler(planarAverageProbe);

    para->setHasWallModelMonitor(true);
    auto wallModelProbe =
        std::make_shared<WallModelProbe>(para, cudaMemoryManager, para->getOutputPath(), "wallModelProbe", 0,
                                         tStartTmpAveraging / deltaT, 100, 0, 100, false, true, true, false, useSurfaceLayer);
    para->addSampler(wallModelProbe);

    for (int probeHeight : { 27, 37, 90, 153, 186, 335 }) {
        std::string name = "planeProbe_" + std::to_string(probeHeight);
        SPtr<Probe> planeProbe =
            std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), name, tStartSampling / deltaT,
                                    tOutProbe / deltaT, tStartOutProbe / deltaT, tOutProbe / deltaT, false, false, true);
        planeProbe->addProbePlane(0, 0, static_cast<real>(probeHeight), lengthX, lengthY, deltaX);
        planeProbe->addAllAvailableStatistics();
        para->addSampler(planeProbe);
    }

    auto tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    Simulation sim(para, cudaMemoryManager, gridBuilderFacade->getGridBuilder(), &bcFactory, tmFactory);
    sim.run();
}

int main(int argc, char* argv[])
{
    if (argv == NULL)
        return 1;

    try {
        vf::logging::Logger::initializeLogger();

        auto config = vf::basics::loadConfig(argc, argv, defaultConfigFile);

        run(config);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cout << "Log initialization failed: " << ex.what() << "\n";
    }

    catch (const std::bad_alloc& e) {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    } catch (const std::exception& e) {
        VF_LOG_CRITICAL("exception: {}", e.what());
    } catch (...) {
        VF_LOG_CRITICAL("Unknown exception!");
    }

    return 0;
}