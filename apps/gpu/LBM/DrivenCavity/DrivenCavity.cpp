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
//! \file LidDrivenCavity.cpp
//! \ingroup Applications
//! \author Martin Schoenherr, Stephan Lenz
//=======================================================================================
#define _USE_MATH_DEFINES
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "Core/LbmOrGks.h"
#include "Core/Logger/Logger.h"
#include "Core/VectorTypes.h"
#include "PointerDefinitions.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

//////////////////////////////////////////////////////////////////////////

// #include "GksMeshAdapter/GksMeshAdapter.h"
// #include "GksGpu/DataBase/DataBase.h"
// #include "GksGpu/Initializer/Initializer.h"
// #include "GksGpu/Parameters/Parameters.h"
// #include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"
// #include "GksGpu/BoundaryConditions/BoundaryCondition.h"
// #include "GksGpu/BoundaryConditions/IsothermalWall.h"
// #include "GksGpu/TimeStepping/NestedTimeStep.h"
// #include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
// #include "GksGpu/Analyzer/CupsAnalyzer.h"
// #include "GksGpu/CudaUtility/CudaUtility.h"
// #include "GksGpu/Output/VtkWriter.h"

//////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    try {
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////
        std::string path("./output/DrivenCavity");
        std::string simulationName("LidDrivenCavity");

        const real L = 1.0;
        const real Re = 1000.0;
        const real velocity = 1.0;
        const real dt = (real)0.5e-3;
        const uint nx = 64;

        const uint timeStepOut = 1000;
        const uint timeStepEnd = 10000;

        // switch between LBM and GKS solver here
        // LbmOrGks lbmOrGks = GKS;
        LbmOrGks lbmOrGks = LBM;

        //////////////////////////////////////////////////////////////////////////
        // setup logger
        //////////////////////////////////////////////////////////////////////////

        logging::Logger::addStream(&std::cout);
        logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
        logging::Logger::timeStamp(logging::Logger::ENABLE);
        logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

        //////////////////////////////////////////////////////////////////////////
        // setup gridGenerator
        //////////////////////////////////////////////////////////////////////////

        auto gridFactory = GridFactory::make();
        gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
        auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////

        real dx = L / real(nx);

        gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L, 0.5 * L, 0.5 * L, 0.5 * L, dx);

        gridBuilder->setPeriodicBoundaryCondition(false, false, false);

        gridBuilder->buildGrids(lbmOrGks, false);

        //////////////////////////////////////////////////////////////////////////
        // branch between LBM and GKS
        //////////////////////////////////////////////////////////////////////////

        if (lbmOrGks == LBM) {
            SPtr<Parameter> para = std::make_shared<Parameter>();
            BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
            vf::gpu::Communicator &communicator = vf::gpu::Communicator::getInstance();

            //////////////////////////////////////////////////////////////////////////
            // compute parameters in lattice units
            //////////////////////////////////////////////////////////////////////////

            const real velocityLB = velocity * dt / dx; // LB units

            const real vxLB = velocityLB / sqrt(2.0); // LB units
            const real vyLB = velocityLB / sqrt(2.0); // LB units

            const real viscosityLB = nx * velocityLB / Re; // LB units

            *logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
            *logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

            //////////////////////////////////////////////////////////////////////////
            // set parameters
            //////////////////////////////////////////////////////////////////////////

            para->setOutputPath(path);
            para->setOutputPrefix(simulationName);

            para->setPrintFiles(true);

            para->setVelocityLB(velocityLB);
            para->setViscosityLB(viscosityLB);

            para->setVelocityRatio(velocity / velocityLB);
            para->setDensityRatio(1.0);

            para->setTimestepOut(timeStepOut);
            para->setTimestepEnd(timeStepEnd);

            para->setMainKernel("CumulantK17CompChimRedesigned");

            //////////////////////////////////////////////////////////////////////////
            // set boundary conditions
            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setNoSlipBoundaryCondition(SideType::PX);
            gridBuilder->setNoSlipBoundaryCondition(SideType::MX);
            gridBuilder->setNoSlipBoundaryCondition(SideType::PY);
            gridBuilder->setNoSlipBoundaryCondition(SideType::MY);
            gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);
            gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);

            bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack);
            bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocitySimpleBounceBackCompressible);

            //////////////////////////////////////////////////////////////////////////
            // set copy mesh to simulation
            //////////////////////////////////////////////////////////////////////////

            auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
            SPtr<GridProvider> gridGenerator =
                GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

            para->setCalcParticles(false);
            para->setCalcMedian(false);
            para->setCalcTurbulenceIntensity(false);
            para->setDoCheckPoint(false);
            para->setUseMeasurePoints(false);
            para->setDiffOn(false);
            para->setCalcPlaneConc(false);
            para->setUseWale(false);
            para->setSimulatePorousMedia(false);
            //////////////////////////////////////////////////////////////////////////
            // run simulation
            //////////////////////////////////////////////////////////////////////////

            Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
            sim.run();
        } // else {
        //     CudaUtility::setCudaDevice(0);

        //     Parameters parameters;

        //     //////////////////////////////////////////////////////////////////////////
        //     // compute remaining parameters
        //     //////////////////////////////////////////////////////////////////////////

        //     const real vx = velocity / sqrt(2.0);
        //     const real vy = velocity / sqrt(2.0);

        //     parameters.K  = 2.0;
        //     parameters.Pr = 1.0;

        //     const real Ma = (real)0.1;

        //     real rho = 1.0;

        //     real cs     = velocity / Ma;
        //     real lambda = c1o2 * ((parameters.K + 5.0) / (parameters.K + 3.0)) / (cs * cs);

        //     const real mu = velocity * L * rho / Re;

        //     *logging::out << logging::Logger::INFO_HIGH << "mu  = " << mu << " m^2/s\n";

        //     *logging::out << logging::Logger::INFO_HIGH << "CFL = " << dt * (velocity + cs) / dx << "\n";

        //     //////////////////////////////////////////////////////////////////////////
        //     // set parameters
        //     //////////////////////////////////////////////////////////////////////////

        //     parameters.mu = mu;

        //     parameters.dt = dt;
        //     parameters.dx = dx;

        //     parameters.lambdaRef = lambda;

        //     //////////////////////////////////////////////////////////////////////////
        //     // set copy mesh to simulation
        //     //////////////////////////////////////////////////////////////////////////

        //     GksMeshAdapter meshAdapter(gridBuilder);

        //     meshAdapter.inputGrid();

        //     auto dataBase = std::make_shared<DataBase>("GPU");

        //     //////////////////////////////////////////////////////////////////////////
        //     // set boundary conditions
        //     //////////////////////////////////////////////////////////////////////////

        //     SPtr<BoundaryCondition> bcLid =
        //         std::make_shared<IsothermalWall>(dataBase, Vec3(vx, vy, 0.0), lambda, false);
        //     SPtr<BoundaryCondition> bcWall =
        //         std::make_shared<IsothermalWall>(dataBase, Vec3(0.0, 0.0, 0.0), lambda, false);

        //     bcLid->findBoundaryCells(meshAdapter, false, [&](Vec3 center) {
        //         return center.z > 0.5 && center.x > -0.5 && center.x < 0.5 && center.y > -0.5 && center.y < 0.5;
        //     });

        //     bcWall->findBoundaryCells(meshAdapter, true, [&](Vec3 center) {
        //         return center.x < -0.5 || center.x > 0.5 || center.y < -0.5 || center.y > 0.5 || center.z < -0.5;
        //     });

        //     dataBase->boundaryConditions.push_back(bcLid);
        //     dataBase->boundaryConditions.push_back(bcWall);

        //     //////////////////////////////////////////////////////////////////////////
        //     // set initial condition and upload mesh and initial condition to GPGPU
        //     //////////////////////////////////////////////////////////////////////////

        //     dataBase->setMesh(meshAdapter);

        //     Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> ConservedVariables {
        //         return toConservedVariables(PrimitiveVariables(rho, 0.0, 0.0, 0.0, lambda), parameters.K);
        //     });

        //     dataBase->copyDataHostToDevice();

        //     Initializer::initializeDataUpdate(dataBase);

        //     VtkWriter::write(dataBase, parameters, path + "/" + simulationName + "_0");

        //     //////////////////////////////////////////////////////////////////////////
        //     // set analyzers
        //     //////////////////////////////////////////////////////////////////////////

        //     CupsAnalyzer cupsAnalyzer(dataBase, false, 60.0, true, 10000);

        //     ConvergenceAnalyzer convergenceAnalyzer(dataBase, 10000);

        //     cupsAnalyzer.start();

        //     //////////////////////////////////////////////////////////////////////////
        //     // run simulation
        //     //////////////////////////////////////////////////////////////////////////

        //     for (uint iter = 1; iter <= timeStepEnd; iter++) {
        //         TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        //         if (iter % timeStepOut == 0) {
        //             dataBase->copyDataDeviceToHost();

        //             VtkWriter::write(dataBase, parameters, path + "/" + simulationName + "_" + std::to_string(iter));
        //         }

        //         int crashCellIndex = dataBase->getCrashCellIndex();
        //         if (crashCellIndex >= 0) {
        //             *logging::out << logging::Logger::LOGGER_ERROR
        //                           << "Simulation crashed at CellIndex = " << crashCellIndex << "\n";
        //             dataBase->copyDataDeviceToHost();
        //             VtkWriter::write(dataBase, parameters, path + "/" + simulationName + "_" + std::to_string(iter));

        //             break;
        //         }

        //         dataBase->getCrashCellIndex();

        //         cupsAnalyzer.run(iter, parameters.dt);

        //         convergenceAnalyzer.run(iter);
        //     }
        // }
    } catch (const std::bad_alloc &e) {

        *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
    } catch (const std::exception &e) {

        *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
    } catch (std::string &s) {

        *logging::out << logging::Logger::LOGGER_ERROR << s << "\n";
    } catch (...) {
        *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
    }

    return 0;
}
