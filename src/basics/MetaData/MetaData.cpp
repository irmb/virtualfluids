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
//! \addtogroup MetaData
//! \ingroup basics
//! \{
//! \author Soeren Peters
//=======================================================================================
#include "MetaData.h"

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <logger/Logger.h>

#include "buildInfo.h"

namespace vf::basics
{

MetaData::MetaData()
{
    simulation.startDateTime = vf::basics::getCurrentTime();

    buildInfo.git_commit_hash = buildInfo::gitCommitHash();
    buildInfo.git_branch = buildInfo::gitBranch();
    buildInfo.build_type = buildInfo::buildType();
    buildInfo.compiler_flags = buildInfo::compilerFlags();
    buildInfo.remote = buildInfo::buildMachine();
    buildInfo.precision = buildInfo::precision();
    buildInfo.compiler_definitions = buildInfo::compilerDefinitions();
    buildInfo.compiler = buildInfo::compiler();
    buildInfo.compiler_version = buildInfo::compiler_version();
#ifdef VF_MPI
    buildInfo.mpi_library = buildInfo::mpi_library();
    buildInfo.mpi_version = buildInfo::mpi_version();
#endif
#ifdef VF_OPENMP
    buildInfo.openmp_library = buildInfo::openmp_library();
    buildInfo.openmp_version = buildInfo::openmp_version();
#endif
}

std::string getCurrentTime()
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    return oss.str();
}

void logPreSimulation(const MetaData& meta_data)
{
    printf("\n");
    VF_LOG_INFO("Start Running {} simulation...", meta_data.name);
    VF_LOG_INFO("Simulation Start Time: {}", meta_data.simulation.startDateTime);
    printf("\n");
    VF_LOG_INFO("world parameter:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("dt [s]                 = {}", meta_data.discretization.dt);
    VF_LOG_INFO("world_length   [m]     = {}", meta_data.world.length);
    VF_LOG_INFO("world_velocity [m/s]   = {}", meta_data.world.velocity);
    VF_LOG_INFO("dx [m]                 = {}", meta_data.discretization.dx);
    printf("\n");
    VF_LOG_INFO("LB parameter:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("Re                     = {}", meta_data.simulation.reynoldsNumber);
    VF_LOG_INFO("lb_velocity [dx/dt]    = {}", meta_data.simulation.lb_velocity);
    VF_LOG_INFO("lb_viscosity [dx^2/dt] = {}", meta_data.simulation.lb_viscosity);
    printf("\n");
    VF_LOG_INFO("simulation parameter:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("number of nodes                = {}", meta_data.discretization.totalNumberOfNodes);
    VF_LOG_INFO("number of level                = {}", meta_data.discretization.numberOfLevels);
    for (uint i = 0; i < meta_data.discretization.numberOfLevels; ++i)
        VF_LOG_INFO("number of nodes on level {} = {}", i, meta_data.discretization.numberOfNodesPerLevel[i]);
    VF_LOG_INFO("total number of time steps     = {}", meta_data.simulation.numberOfTimeSteps);
    VF_LOG_INFO("collision kernel               = {}", meta_data.simulation.collisionKernel);
    VF_LOG_INFO("quadric limiter                = {}, {}, {}", meta_data.simulation.quadricLimiters[0],
                meta_data.simulation.quadricLimiters[1], meta_data.simulation.quadricLimiters[2]);
    printf("\n");
    VF_LOG_INFO("Build Info:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("git commit hash: {}", meta_data.buildInfo.git_commit_hash);
    VF_LOG_INFO("git branch: {}", meta_data.buildInfo.git_branch);
    VF_LOG_INFO("build type: {}", meta_data.buildInfo.build_type);
    VF_LOG_INFO("remote: {}", meta_data.buildInfo.remote);
    VF_LOG_INFO("Precision: {}", meta_data.buildInfo.precision);
    VF_LOG_INFO("compiler: {}", meta_data.buildInfo.compiler);
    VF_LOG_INFO("compiler version: {}", meta_data.buildInfo.compiler_version);
    VF_LOG_INFO("compiler flags: {}", meta_data.buildInfo.compiler_flags);
    VF_LOG_INFO("compiler definitions: {}", meta_data.buildInfo.compiler_definitions);
    printf("\n");
    VF_LOG_INFO("Hardware Info:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("Number of Processes: {}", meta_data.numberOfProcesses);
    VF_LOG_INFO("Number of Threads: {}", meta_data.numberOfThreads);
    VF_LOG_INFO("hardware used: {}", meta_data.vf_hardware);
    printf("\n");
    if (meta_data.vf_hardware == "GPU") {
        VF_LOG_INFO("GPU Info:");
        VF_LOG_INFO("--------------");
        for (const auto& gpu : meta_data.gpus) {
            VF_LOG_INFO("name: {}, compute capability: {}", gpu.name, gpu.compute_capability);
        }
    }
    printf("\n");
}

void logPostSimulation(const MetaData& meta_data)
{
    printf("\n");
    VF_LOG_INFO("... finish Running simulation...");
    VF_LOG_INFO("Total runtime: {:.0f} ms", meta_data.simulation.runtimeSeconds * 1000);
    VF_LOG_INFO("NUPS: {:.0f} ({:04.1f}e6 NUPS)", meta_data.simulation.nups, meta_data.simulation.nups * 1e-6);
    printf("\n");
}

} // namespace vf::basics

//! \}
