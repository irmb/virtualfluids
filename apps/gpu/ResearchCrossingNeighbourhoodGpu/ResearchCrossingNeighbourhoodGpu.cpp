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
//! \author Hussein Alihussein
//! \brief GPU app for ResearchCrossingNeighbourhood grid generation and winding diagnostics.
//! \note Generated with the assistance of OpenAI Codex (GPT-5.1 Codex Max). Reviewed and adapted by author.
//=======================================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>
#include <basics/Timer/Timer.h>
#include <basics/geometry3d/GbCuboid3D.h>
#include <basics/geometry3d/GbSystem3D.h>
#include <basics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/geometry3d/winding/GridWindingFastDefaults.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

#include <filesystem>
#include <optional>

#include "GridGenerator/geometries/BoundingBox/BoundingBox.h"
#include "GridGenerator/io/GridWindingDiagnosticsGpu.h"

#include "GridGenerator/geometries/Triangle/Triangle.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/geometries/Vertex/Vertex.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridImp.h"
#include "GridGenerator/grid/GridDimensions.h"
#include "GridGenerator/grid/MultipleGridBuilderFacade.h"

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Parameter/Parameter.h"

#include "utilities/transformator/TransformatorImp.h"

#include <basics/geometry3d/Axis.h>

using namespace vf::basics::constant;
using namespace vf::gpu;
namespace
{
using vf::basics::ConfigurationFile;

struct GridLayout
{
    std::array<real, 3> domainMin{};
    std::array<real, 3> domainMax{};
    std::array<real, 3> ghostMin{};
    std::array<real, 3> ghostMax{};
    std::array<uint, 3> nodeCounts{};
};

GridLayout computeGridLayout(const std::array<real, 3>& domainMin,
                             const std::array<real, 3>& domainMax,
                             real deltaX)
{
    GridLayout layout{};
    layout.domainMin = domainMin;
    layout.domainMax = domainMax;

    const real halfDelta = static_cast<real>(0.5) * deltaX;
    for (int axis = 0; axis < 3; ++axis) {
        const real min = layout.domainMin[axis];
        const real max = layout.domainMax[axis];
        if (max <= min) {
            std::ostringstream message;
            message << "Invalid domain extent along axis " << axis << ": min=" << min << ", max=" << max;
            throw std::invalid_argument(message.str());
        }

        layout.ghostMin[axis] = min - halfDelta;
        layout.ghostMax[axis] = max + halfDelta;

        const real extent = max - min;
        const double cells = static_cast<double>(extent / deltaX);
        const long rounded = std::lround(cells);
        if (rounded < 1)
            throw std::invalid_argument("Domain extent too small for requested deltaX.");
        layout.nodeCounts[axis] = static_cast<uint>(rounded) + 2u;
    }

    return layout;
}

std::shared_ptr<TriangularMesh> loadAndTransformTriMesh(const std::string& filename, real shiftX, real shiftY, real shiftZ)
{
    auto mesh = std::make_shared<TriangularMesh>(filename);

    TransformatorImp transform(1.0, Vertex(shiftX, shiftY, shiftZ));
    transform.transformWorldToGrid(*mesh);

    BoundingBox bbox = BoundingBox::makeInvalidMinMaxBox();
    for (const Triangle& triangle : mesh->triangleVec) {
        bbox.setMinMax(triangle);
    }
    mesh->setMinMax(bbox);

    return mesh;
}

void run(const ConfigurationFile& config)
{
    const std::string path = config.getValue<std::string>("path");
    const bool writeMissingLinks = config.getValue<bool>("WriteMissingLinks");
    const bool writeQLines = config.getValue<bool>("WriteQLines");

    std::error_code outputDirStatus;
    std::filesystem::create_directories(path, outputDirStatus);
    (void)outputDirStatus;

    const std::string cityFilename = config.getValue<std::string>("CityFilename");
    const std::string bodenFilename = config.getValue<std::string>("BodenFilename");

    const real Re = config.getValue<real>("Re"); 
    
    const real velocity = config.getValue<real>("velocity"); 
    const real deltaX = config.getValue<real>("deltaX");
#if defined(VF_HAS_FAST_WINDING)
    const real fastWindingThreshold = config.getValue<real>("FastWindingThreshold");
#endif
    const uint timeStepOut = config.getValue<uint>("timeStepOut");
    const uint timeStepEnd = config.getValue<uint>("timeStepEnd");

    auto city = std::make_shared<GbTriFaceMesh3D>();
    auto bodenSurface = std::make_shared<GbTriFaceMesh3D>();
    bool removeRedundant = false;
    city->readMeshFromSTLFileBinary(cityFilename, removeRedundant);
    bodenSurface->readMeshFromSTLFileBinary(bodenFilename, removeRedundant);

    const real initialX = city->getX1Minimum();
    const real initialY = city->getX2Minimum();
    const real initialZ = city->getX3Minimum();
    city->translate(-initialX, -initialY, -initialZ);
    bodenSurface->translate(-initialX, -initialY, -initialZ);

    gb_system_3d::writeGeoObject(city.get(), path + "/geo/city", WbWriterVtkXmlBinary::getInstance());
    gb_system_3d::writeGeoObject(bodenSurface.get(), path + "/geo/boden", WbWriterVtkXmlBinary::getInstance());

    const std::vector<real> gridMinValues = config.getVector<real>("GridCubeMin");
    const std::vector<real> gridMaxValues = config.getVector<real>("GridCubeMax");
    if (gridMinValues.size() != 3 || gridMaxValues.size() != 3) {
        throw std::runtime_error("GridCubeMin and GridCubeMax must have 3 entries each.");
    }
    const std::array<real, 3> gridMin{ gridMinValues[0], gridMinValues[1], gridMinValues[2] };
    const std::array<real, 3> gridMax{ gridMaxValues[0], gridMaxValues[1], gridMaxValues[2] };
    const auto layout = computeGridLayout(gridMin, gridMax, deltaX);

    GbCuboid3DPtr gridCube(new GbCuboid3D(layout.domainMin[0], layout.domainMin[1], layout.domainMin[2],
                                          layout.domainMax[0], layout.domainMax[1], layout.domainMax[2]));
    gb_system_3d::writeGeoObject(gridCube.get(), path + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

    const real velocityLB = velocity;
    const std::array<real, 3> domainSize{
        layout.domainMax[0] - layout.domainMin[0],
        layout.domainMax[1] - layout.domainMin[1],
        layout.domainMax[2] - layout.domainMin[2]
    };
    const real nx = static_cast<real>(layout.nodeCounts[0] - 2u);
    const real viscosityLB = nx * velocityLB / Re;

    VF_LOG_INFO("ResearchCrossingNeighbourhoodGpu:");
    VF_LOG_INFO("Domain extent (m) = {} x {} x {}", domainSize[0], domainSize[1], domainSize[2]);
    VF_LOG_INFO("Domain cells = {} x {} x {}", static_cast<unsigned long long>(layout.nodeCounts[0] - 2u),
                static_cast<unsigned long long>(layout.nodeCounts[1] - 2u),
                static_cast<unsigned long long>(layout.nodeCounts[2] - 2u));
    VF_LOG_INFO("velocity    = {}", velocity);
    VF_LOG_INFO("velocityLB  = {}", velocityLB);
    VF_LOG_INFO("viscosityLB = {}", viscosityLB);
    VF_LOG_INFO("Re = {}", Re);
    VF_LOG_INFO("dx = {}", deltaX);
    VF_LOG_INFO("Ghosted grid min=({}, {}, {}), max=({}, {}, {}), nodes=({}, {}, {})", layout.ghostMin[0],
                layout.ghostMin[1], layout.ghostMin[2], layout.ghostMax[0], layout.ghostMax[1], layout.ghostMax[2],
                layout.nodeCounts[0], layout.nodeCounts[1], layout.nodeCounts[2]);

    VF_LOG_TRACE("Preparing city mesh for grid generation.");
    auto cityMesh = loadAndTransformTriMesh(cityFilename, -initialX, -initialY, -initialZ);
    VF_LOG_TRACE("City mesh loaded, generating GbTriFaceMesh3D representation.");
    cityMesh->generateGbTriFaceMesh3D();
    VF_LOG_TRACE("City mesh GbTriFaceMesh3D generated.");
    VF_LOG_TRACE("Preparing Boden mesh for grid generation.");
    auto bodenMesh = loadAndTransformTriMesh(bodenFilename, -initialX, -initialY, -initialZ);
    VF_LOG_TRACE("Boden mesh loaded, generating GbTriFaceMesh3D representation.");
    bodenMesh->generateGbTriFaceMesh3D();
    VF_LOG_TRACE("Boden mesh GbTriFaceMesh3D generated.");

    auto comm = vf::parallel::MPICommunicator::getInstance();
    vf::parallel::Communicator& communicator = *comm;
    const bool isRoot = communicator.isRoot();
    const int numberOfProcesses = communicator.getNumberOfProcesses();
    const int processId = communicator.getProcessID();

    const bool useMultiGpuConfig = config.getValue<bool>("UseMultiGPU");
    bool multiGpuAvailable = numberOfProcesses > 1;
    bool useMultiGpu = useMultiGpuConfig && multiGpuAvailable;
    if (!useMultiGpuConfig && isRoot)
        VF_LOG_INFO("UseMultiGPU disabled via configuration file.");
    if (multiGpuAvailable && !useMultiGpu && isRoot) {
        VF_LOG_WARNING("MPI launch detected but multi-GPU domain decomposition disabled; each rank will process the full domain.");
    }
    if (!multiGpuAvailable && useMultiGpuConfig && isRoot) {
        VF_LOG_INFO("UseMultiGPU requested but running on a single rank; falling back to single-GPU mode.");
    }
    if (isRoot) {
        VF_LOG_INFO("Multi-GPU mode {}", (useMultiGpu ? "enabled" : "disabled"));
    }

    auto domainDimensions = std::make_shared<GridDimensions>();
    domainDimensions->minX = layout.domainMin[0];
    domainDimensions->maxX = layout.domainMax[0];
    domainDimensions->minY = layout.domainMin[1];
    domainDimensions->maxY = layout.domainMax[1];
    domainDimensions->minZ = layout.domainMin[2];
    domainDimensions->maxZ = layout.domainMax[2];
    domainDimensions->delta = deltaX;

    auto gridBuilderCore = std::make_shared<MultipleGridBuilder>();
    // MultipleGridBuilder defaults to POINT_IN_OBJECT. For this app, prefer FAST_WINDING when compiled in.
    // Other available strategies in GridFactory are RAYCASTING and POINT_UNDER_TRIANGLE.
#if defined(VF_HAS_FAST_WINDING)
    // Inject a configured strategy instance to keep MultipleGridBuilder free of fast-winding-specific settings.
    gridBuilderCore->setTriangularMeshDiscretizationStrategy(
        std::make_shared<FastWindingDiscretizationStrategy>(
            vf::grid_winding::FastWindingDefaultAccuracyScale,
            static_cast<float>(fastWindingThreshold),
            vf::grid_winding::FastWindingDefaultTolerance));
#endif
    auto gridBuilderFacade =
        std::make_unique<MultipleGridBuilderFacade>(gridBuilderCore, domainDimensions, std::optional<real>(deltaX));

    VF_LOG_TRACE("Adding geometry to builder facade.");
    gridBuilderFacade->addGeometry(cityMesh);
    gridBuilderFacade->addGeometry(bodenMesh);
    gridBuilderFacade->setPeriodicBoundaryCondition(true, true, false);

    if (useMultiGpu && numberOfProcesses > 1) {
        const real spanX = layout.domainMax[0] - layout.domainMin[0];
        const real minSplit = layout.domainMin[0];
        if (spanX <= static_cast<real>(0.0)) {
            VF_LOG_WARNING("Domain span in X is zero; falling back to single-GPU grid.");
        } else {
            for (int split = 1; split < numberOfProcesses; ++split) {
                const real coord =
                    minSplit + spanX * (static_cast<real>(split) / static_cast<real>(numberOfProcesses));
                gridBuilderFacade->addDomainSplit(coord, Axis::x);
            }
        }
    }

    const uint subdomainIndex = useMultiGpu ? static_cast<uint>(processId) : 0u;
    VF_LOG_TRACE("Creating grids for subdomain {}.", subdomainIndex);
    gridBuilderFacade->createGrids(subdomainIndex);

    VF_LOG_TRACE("Applying boundary conditions via facade.");
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::PZ, velocityLB, 0.2 * velocityLB, 0.0);
    gridBuilderFacade->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

    auto gridBuilder = gridBuilderFacade->getGridBuilder();
    const auto grids = gridBuilder->getGrids();
    GridScalingFactory scalingFactory;
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);
    vf::gpu::grid_winding::runBoundaryDiagnosticsForSurfaces(grids,
                                                             city,
                                                             bodenSurface,
                                                             writeMissingLinks,
                                                             writeQLines,
                                                             path,
                                                             comm,
                                                             isRoot);

    auto para = std::make_shared<Parameter>(numberOfProcesses, processId);
    para->setOutputPath(path);
    para->setOutputPrefix("ResearchCrossingNeighbourhoodGpu");
    para->setPrintFiles(true);
    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(1.0);
    para->setDensityRatio(1.0);
    para->setTimestepOut(timeStepOut);
    para->setTimestepEnd(timeStepEnd);
    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);
    para->setInitialCondition([](real /*x*/, real /*y*/, real /*z*/, real& rho, real& vx, real& vy, real& vz) {
        rho = c0o1;
        vx  = static_cast<real>(0.1);
        vy  = static_cast<real>(0.1);
        vz  = static_cast<real>(0.0);
    });
    para->setAllNodesAllFeatures(false);
    para->useReducedCommunicationAfterFtoC = useMultiGpu;

    BoundaryConditionFactory bcFactory;
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    VF_LOG_INFO("Starting ResearchCrossingNeighbourhoodGpu simulation.");

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);
    sim.run();

    VF_LOG_INFO("Finished ResearchCrossingNeighbourhoodGpu simulation.");
}
} // namespace

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();

        const std::string defaultConfigPath = "./apps/gpu/ResearchCrossingNeighbourhoodGpu/ResearchCrossingNeighbourhoodGpu.cfg";
        std::string configPath = defaultConfigPath;
        bool configExplicit = false;

        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);
            if (arg == "--config") {
                if (i + 1 >= argc)
                    throw std::invalid_argument("--config requires a value");
                configPath = argv[++i];
                configExplicit = true;
            } else if (arg.rfind("--", 0) == 0) {
                VF_LOG_WARNING("Unrecognized option {} ignored", arg);
            } else if (!configExplicit) {
                configPath = arg;
                configExplicit = true;
            } else {
                VF_LOG_WARNING("Ignoring extra positional argument {}", arg);
            }
        }

        ConfigurationFile config;
        config.load(configPath);
        run(config);
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}
