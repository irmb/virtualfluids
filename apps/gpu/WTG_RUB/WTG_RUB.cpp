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
//! \file WTG_RUB.cpp
//! \ingroup Applications
//! \author Martin Schoenherr
//=======================================================================================
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>
#include <filesystem>

//////////////////////////////////////////////////////////////////////////

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/GridFactory.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"

#include <parallel/MPICommunicator.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// const real L  = 1.0;

const real velocity  = 1.0;

int variant = 1;
real rotationOfCity;
real z_offset = 0.0; // only if baseplate is in use (currently not!! not important)
int dataN;
std::vector<real> dataZ;
std::vector<real> dataVelocity;
std::string simulationName("");

// 1: original setup of Lennard Lux (6 level, 4.0 cm -> 1.25 mm)
// 2: setup 1 of MSch               (4 level, 1.0 cm -> 1.25 mm)
// 3: setup 2 of MSch               (5 level, 1.6 cm -> 1.0  mm)
// 4: setup 3 of MSch (small/test)  (3 level, 4.0 cm -> 1.0  cm)
int setupDomain = 4;

// std::string path("D:/out/WTG_RUB"); // Mollok
// std::string inputPath("D:/out/WTG_RUB/input/");

std::string path("/workspaces/VirtualFluids_dev/output/WTG_RUB_Results/"); // Aragorn
std::string inputPath("/workspaces/VirtualFluids_dev/stl/WTG_RUB/");

// const uint timeStepStartOut = 0;
const uint timeStepOut = 10000;
const uint timeStepEnd = 100000;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void convertMPFile(SPtr <MultipleGridBuilder> gridBuilder, real& phi, std::vector<real>& origin, bool& measureVeloProfilesOnly, uint& maxLevel);

void addFineGrids(SPtr<MultipleGridBuilder> gridBuilder, uint &maxLevel, real &rotationOfCity);

void readVelocityProfile();

std::string chooseVariation();

void multipleLevel(const std::string& configPath)
{
    vf::parallel::Communicator &communicator = *vf::parallel::MPICommunicator::getInstance();
    //auto gridBuilder = MultipleGridBuilder::makeShared();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real dx = 0;
    real viscosityLB = (real)1.0e-03;
    uint maxLevel    = 1;

    if (setupDomain == 1) {
        dx          = (real)4;
        maxLevel    = 5;
        viscosityLB = (real)3.75e-06; // LB units
    } else if (setupDomain == 2) {
        dx          = (real)1;
        maxLevel    = 3;
        viscosityLB = (real)1.5e-05; // LB units
    } else if (setupDomain == 3) {
        dx          = (real)1.6;
        maxLevel    = 4;
        viscosityLB = (real)9.375e-06; // LB units
    } else if (setupDomain == 4) {
        dx = (real)4.0;
        maxLevel = 2;
        viscosityLB = (real)3.75e-06; // LB units
    }

    real x_min = 0.0;
    real x_max = 1250.0;
    real y_min = 0.0;
    real y_max = 190.0;
    real z_min = 0.0 + z_offset;
    real z_max = 160.0 + z_offset;

    //auto RubSTL      = std::make_shared<TriangularMesh>(inputPath + "stl/Var02_0deg_FD_b.stl");
    auto RubSTL      = std::make_shared<TriangularMesh>(inputPath + "stl/" + chooseVariation() + ".stl");
    std::vector<real> originOfCityXY = { 600.0, y_max / 2, z_offset };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // OPTIONS

    int maxLevelAtInflow = maxLevel;
    // there are max. 3 levels (0,1,2) at inflow-/x_min-Side if maxLevel >= 3; OTHERWISE: maxLevel
    if (maxLevel >= 3)
        maxLevelAtInflow = 2;

    // MeasurePoints [MP01-15: lvl maxLevel],[MP16-41: lvl 1]; disable when reducing numberOfLevels --> dx might be too
    // large if MP01-15 are used with low resolution dx, MPs might be placed in solid City-geometry
    bool useMP                   = false;//true;
    bool measureVeloProfilesOnly = false;

    // Two Components: true->DiffOn, false->DiffOff
    // bool diffOnOff = false;

    // Resetting diff or flow field, e.g. after restart, do not reset diff/flow at start of measureRun ;-)
    // bool reset_diff = false;
    bool reset_flow = false;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    gridBuilder->addCoarseGrid(x_min, y_min, z_min,
                               x_max, y_max, z_max, dx);

    gridBuilder->setNumberOfLayers(0, 0);

    addFineGrids(gridBuilder, maxLevel, rotationOfCity);

    //// adding solid CityGeometry to gridbuilder
    gridBuilder->addGeometry(RubSTL);

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(false); // buildGrids() has to be called before setting the BCs!!!!

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vf::basics::ConfigurationFile config;
    config.load(configPath);

    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcesses(), communicator.getProcessID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const real velocityLB = (real)0.0844; // LB units

    //const real vx = velocityLB / (real)sqrt(2.0); // LB units
    //const real vy = velocityLB / (real)sqrt(2.0); // LB units

    VF_LOG_INFO("velocityLB [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("viscosityLB [dx^2/dt] = {}", viscosityLB);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setDevices(std::vector<uint>{(uint)0});

    para->setOutputPath( path );
    para->setOutputPrefix( "Unified_" + simulationName );

    para->setPrintFiles(true);

    para->setMaxLevel(maxLevel);

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);

    para->setVelocityRatio(velocity/ velocityLB);

    para->configureMainKernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = (real)0.0; //(6 * velocityLB * coordZ * (L - coordZ) / (L * L));
        vy  = (real)0.0;
        vz  = (real)0.0;
    });

    para->setTimestepOut( timeStepOut );
    para->setTimestepEnd( timeStepEnd );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityWithPressureInterpolatedCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    readVelocityProfile();
    for (int level = 0; level <= maxLevelAtInflow; level++) {
        auto inflowBC = std::dynamic_pointer_cast<VelocityBoundaryCondition>(
            gridBuilder->getBoundaryCondition(SideType::MX, level));

        // lambda function [capture list](parameter list){function body}
        // for every z, loop is being run;
        if (inflowBC) {
            inflowBC->setVelocityProfile(
                gridBuilder->getGrid(level), [&](real x, real y, real z, real &vx, real &vy, real &vz) {
                    int i;
                    for (i = 0; i < dataN; i++) {
                        if ((z < dataZ[i]) || (z == dataZ[i]))
                            break;
                    } // remembers i

                    // if z is below current data point --> interpolation between previous and current datapoint
                    if (z < dataZ[i]) {
                        vx = velocityLB * (dataVelocity[i] + (dataVelocity[i + 1] - dataVelocity[i]) /
                                                                 (dataZ[i + 1] - dataZ[i]) * (z - dataZ[i]));
                    } else if (z == dataZ[i]) {
                        vx = velocityLB * dataVelocity[i];
                    }

                    // vx << std::endl; vx = velocityLB;
                    vy = 0.0;
                    vz = 0.0;
                });
        }
    }

    // Resetting Velocity Profile (returning to intial state; do not combine with restart)
    if (reset_flow) {
        para->setInitialCondition([&](real x, real y, real z, real &rho, real &vx, real &vy, real &vz) {
            int i;
            for (i = 0; i < dataN; i++) {
                if ((z < dataZ[i]) || (z == dataZ[i]))
                    break;
            }
            if (z < dataZ[i]) {
                vx = velocityLB * 0.588 *
                     (dataVelocity[i] +
                      (dataVelocity[i + 1] - dataVelocity[i]) / (dataZ[i + 1] - dataZ[i]) * (z - dataZ[i]));
            } else if (z == dataZ[i]) {
                vx = velocityLB * 0.588 * dataVelocity[i];
            }

            // std::cout << "LINE : " << __LINE__ << "\tdataZ = " << dataZ[i] << "\tz = " << z << "\tvx = " << vx <<
            // std::endl;

            // vx = velocityLB;
            vy  = 0.0;
            vz  = 0.0;
            rho = 0.0;
        });
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Intializing MeasurePoints
    para->setUseMeasurePoints(useMP);
    if (para->getUseMeasurePoints()) {
        convertMPFile(gridBuilder, rotationOfCity, originOfCityXY, measureVeloProfilesOnly, maxLevel);
        // Number log-Files for each MeasurePoint: numberOfMPFiles = timeStepEnd/ClockCycle
        para->setclockCycleForMP(timeStepEnd);
        // Number of  logged Timesteps for each file
        para->settimestepForMP(timeStepOut / 100);
        // para->settimestepForMP(timeStepOut);
        para->setmeasurePoints(inputPath + "measurePoints.dat");
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    auto gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
    sim.run();
}

void readVelocityProfile()
{
    // reads velocityProfile.txt, containing values for relative velocity (u_h) in relation to height (z). also fills
    // dataZ,dataVelocity vectors

    std::ifstream inFile;
    inFile.open(inputPath + "VeloProfile.txt");
    //inFile.open(inputPath + "velocityProfile.txt");

    if (inFile.fail()) {
        std::cerr << "Error opening File" << std::endl;
        exit(1);
    }

    int z;
    real velocity;

    // read
    inFile >> dataN;
    for (int k = 0; k < dataN; k++) {
        inFile >> z;
        dataZ.push_back(z + z_offset);

        inFile >> velocity;
        dataVelocity.push_back(velocity);
    }
    inFile.close();
}

void addFineGrids(SPtr<MultipleGridBuilder> gridBuilder, uint &maxLevel, real &rotationOfCity)
{
    if (setupDomain == 1) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
        // creates Cuboids (FG1 to FG3, lvl 1 to lvl 3) and add STLs (FG4 to FG5, lvl 4 to lvl 5) depending on maxLevel
        // and rotationOfCity; also adds FineGrids(FGs) to gridbuilder

        // GridList(CG = coarse grid, fg = fine grid)
        // CG  -> dx = 4 cm;      lvl 0
        // FG1 -> dx = 2 cm;      lvl 1
        // FG2 -> dx = 1 cm;      lvl 2
        // FG3 -> dx = 5 mm;      lvl 3
        // FG4 -> dx = 2,5 mm;    lvl 4
        // FG5 -> dx = 1,25 mm;   lvl 5
        //
        // FineGrid Level 1 ->dx = 2 cm; lvl 1
        auto FG1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 800, 200, 75 + z_offset);

        // FineGrid Level 2 -> dx = 1 cm; lvl 2
        auto FG2_1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG2_2 = std::make_shared<Cuboid>(500, -20,  5 + z_offset, 680, 210, 50 + z_offset);
        auto FG2   = std::make_shared<Conglomerate>();
        FG2->add(FG2_1);
        FG2->add(FG2_2);

        // FineGrid Level 3 ->dx = 5 mm; lvl 3
        auto FG3_1 = std::make_shared<Cuboid>(517, -20, -5 + z_offset, 665, 200, 30 + z_offset);
        auto FG3_2 = std::make_shared<Cuboid>(550, 58, -5 + z_offset, 650, 132, 40 + z_offset);
        auto FG3   = std::make_shared<Conglomerate>();
        FG3->add(FG3_1);
        FG3->add(FG3_2);

        // Adding FineGrids 1 to 5 depending on maxLevel, FG4 and FG5 require different STL-files depending on
        // rotationOfCity
        if (maxLevel >= 1) {
            gridBuilder->addGrid(FG1, 1);
            if (maxLevel >= 2) {
                gridBuilder->addGrid(FG2, 2);
                if (maxLevel >= 3) {
                    gridBuilder->addGrid(FG3, 3);
                    if (maxLevel >= 4) {
                        if (rotationOfCity == 0.0) {
                            auto FG4 = std::make_shared<TriangularMesh>(inputPath + "stl/FG4_0deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        } else {
                            auto FG4 = std::make_shared<TriangularMesh>(inputPath + "stl/FG4_63deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        }

                        if (maxLevel == 5) {
                            if (rotationOfCity == 0.0) {
                                auto FG5 = std::make_shared<TriangularMesh>(inputPath + "stl/FG5_0deg.stl");
                                gridBuilder->addGrid(FG5, 5);
                            } else {
                                auto FG5 = std::make_shared<TriangularMesh>(inputPath + "stl/FG5_63deg.stl");
                                gridBuilder->addGrid(FG5, 5);
                            }
                        }
                    }
                }
            }
        }
    }
    else if (setupDomain == 2) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
        // creates Cuboids (FG1, lvl 1) and add STLs (FG2 to FG3, lvl 2 to lvl 3) depending on maxLevel
        // and rotationOfCity; also adds FineGrids(FGs) to gridbuilder
        //
        // GridList(CG = coarse grid, fg = fine grid)
        // CG  -> dx = 1 cm;      lvl 0
        // FG1 -> dx = 5 mm;      lvl 1
        // FG2 -> dx = 2,5 mm;    lvl 2
        // FG3 -> dx = 1,25 mm;   lvl 3
        //
        // FineGrid Level 1 -> dx = 5 mm; lvl 1
        //auto FG1_1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG1_1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 760, 200, 20 + z_offset);
        auto FG1_2 = std::make_shared<Cuboid>(500, -20,  5 + z_offset, 680, 210, 50 + z_offset);
        auto FG1   = std::make_shared<Conglomerate>();
        FG1->add(FG1_1);
        FG1->add(FG1_2);

        // Adding FineGrids 1 to 5 depending on maxLevel, FG4 and FG5 require different STL-files depending on
        // rotationOfCity
        if (maxLevel >= 1) {
            gridBuilder->addGrid(FG1, 1);
            if (maxLevel >= 2) {
                if (rotationOfCity == 0.0) {
                    auto FG2 = std::make_shared<TriangularMesh>(inputPath + "stl/FG4_0deg.stl");
                    gridBuilder->addGrid(FG2, 2);
                } else {
                    auto FG2 = std::make_shared<TriangularMesh>(inputPath + "stl/FG4_63deg.stl");
                    gridBuilder->addGrid(FG2, 2);
                }

                if (maxLevel == 3) {
                    if (rotationOfCity == 0.0) {
                        auto FG3 = std::make_shared<TriangularMesh>(inputPath + "stl/FG5_0deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    } else {
                        auto FG3 = std::make_shared<TriangularMesh>(inputPath + "stl/FG5_63deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    }
                }
            }
        }
    }
    else if (setupDomain == 3) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
        // creates Cuboids (FG1 to FG2, lvl 1 to lvl 2) and add STLs (FG3 to FG4, lvl 3 to lvl 4) depending on maxLevel
        // and rotationOfCity; also adds FineGrids(FGs) to gridbuilder

        // GridList(CG = coarse grid, fg = fine grid)
        // CG  -> dx = 1.6 cm;   lvl 0
        // FG1 -> dx = 8.0 mm;   lvl 1
        // FG2 -> dx = 4.0 mm;   lvl 2
        // FG3 -> dx = 2.0 mm;   lvl 3
        // FG4 -> dx = 1.0 mm;   lvl 4
        //
        //// FineGrid Level 1 ->dx = 8.0 mm; lvl 1
        // auto FG1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 800, 200, 75 + z_offset);

        // FineGrid Level 1 -> dx = 8.0 mm; lvl 1
        auto FG1_1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 780, 200, 30 + z_offset);
        auto FG1_2 = std::make_shared<Cuboid>(500, -20, 5 + z_offset, 720, 210, 75 + z_offset);
        auto FG1   = std::make_shared<Conglomerate>();
        FG1->add(FG1_1);
        FG1->add(FG1_2);

        // FineGrid Level 2 -> dx = 4.0 mm; lvl 2
        auto FG2_1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG2_2 = std::make_shared<Cuboid>(520, -20, 5 + z_offset, 700, 210, 50 + z_offset);
        auto FG2   = std::make_shared<Conglomerate>();
        FG2->add(FG2_1);
        FG2->add(FG2_2);

        // Adding FineGrids 1 to 4 depending on maxLevel, FG3 and FG4 require different STL-files depending on
        // rotationOfCity
        if (maxLevel >= 1) {
            gridBuilder->addGrid(FG1, 1);
            if (maxLevel >= 2) {
                gridBuilder->addGrid(FG2, 2);
                if (maxLevel >= 3) {
                    if (rotationOfCity == 0.0) {
                        auto FG3 = std::make_shared<TriangularMesh>(inputPath + "stl/FG4_0deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    } else {
                        auto FG3 = std::make_shared<TriangularMesh>(inputPath + "stl/FG4_63deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    }

                    if (maxLevel == 4) {
                        if (rotationOfCity == 0.0) {
                            auto FG4 = std::make_shared<TriangularMesh>(inputPath + "stl/FG5_0deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        } else {
                            auto FG4 = std::make_shared<TriangularMesh>(inputPath + "stl/FG5_63deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        }
                    }
                }
            }
        }
    }
    else if (setupDomain == 4) {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
        // creates Cuboids (FG1 to FG2, lvl 1 to lvl 2) depending on
        // maxLevel and rotationOfCity; also adds FineGrids(FGs) to gridbuilder

        // GridList(CG = coarse grid, fg = fine grid)
        // CG  -> dx = 4 cm;      lvl 0
        // FG1 -> dx = 2 cm;      lvl 1
        // FG2 -> dx = 1 cm;      lvl 2
        //
        // FineGrid Level 1 ->dx = 2 cm; lvl 1
        auto FG1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 800, 200, 75 + z_offset);

        // FineGrid Level 2 -> dx = 1 cm; lvl 2
        auto FG2_1 = std::make_shared<Cuboid>(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG2_2 = std::make_shared<Cuboid>(500, -20, 5 + z_offset, 680, 210, 50 + z_offset);
        auto FG2   = std::make_shared<Conglomerate>();
        FG2->add(FG2_1);
        FG2->add(FG2_2);

        // Adding FineGrids 1 to 2 depending on maxLevel
        if (maxLevel >= 1) {
            gridBuilder->addGrid(FG1, 1);
            if (maxLevel >= 2) {
                gridBuilder->addGrid(FG2, 2);
            }
        }
    }



}

void convertMPFile(SPtr<MultipleGridBuilder> gridBuilder, real &phi, std::vector<real> &originXY, bool &measureVeloProfilesOnly, uint &maxLevel)
{
    // File Reader&Writer for converting MP-Coordinates to Index: MeasurePoint placement requires "measurePoints.dat"
    // with name, node-ID and level. This function can read a txt-File containing the name, X-Y-Z-Coordinates and level
    // of measurePoints. After reading the txt-File and converting X-Y-Z to the node-ID, it writes "measurePoints.dat".
    // Justification for this function: Human Readability and no changes in measurepoint core functions

    // File Opening Procedure
    std::ifstream inFile;
    if (measureVeloProfilesOnly)
        inFile.open(inputPath + "measurePoints_veloProfiles.txt");
    else
        inFile.open(inputPath + "measurePoints.txt");

    // Check for error
    if (inFile.fail()) {
        std::cerr << "Error opening File" << std::endl;
        exit(1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reading Procedure
    std::cout << "phi in degrees:" << phi << std::endl;
    phi = phi * M_PI / 180;
    std::cout << "phi in radians:" << phi << std::endl;

    std::vector<std::string> MP_name;
    std::vector<real> MP_X, MP_Y, MP_Z;
    std::vector<int> MP_level, MP_k;

    std::string name;
    real X, Y, Z;
    uint level, numberOfMeasurePoints;

    inFile >> numberOfMeasurePoints;
    std::cout << "numberOfMeasurePoints: " << numberOfMeasurePoints << " ";
    std::cout << "Coordinates from File\n";
    for (uint k = 0; k < numberOfMeasurePoints; k++) {
        inFile >> name;
        MP_name.push_back(name);
        std::cout << "Name: " << MP_name[k] << " ";

        inFile >> X;
        MP_X.push_back(X);
        std::cout << "\t\tX: " << MP_X[k] << " ";

        inFile >> Y;
        MP_Y.push_back(Y);
        std::cout << "\t\tY: " << MP_Y[k] << " ";

        inFile >> Z;
        if (((variant > 3 && variant < 7) || (variant > 9 && variant <= 12)) && k == 14)
            Z += 2.25; // account for angled roof
        MP_Z.push_back(Z);
        std::cout << "\t\tZ: " << MP_Z[k] + z_offset << " ";

        inFile >> level;
        if (level > maxLevel)
            level = maxLevel;
        MP_level.push_back(level);
        std::cout << "\t\tLevel: " << MP_level[k] << std::endl;
    }
    inFile.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real X_temp, Y_temp;
    // Transformation for phi radians around centre of city for MP[1..15]
    if (!phi == 0) {
        std::cout << "Calculating new Coordinates for MP01 to MP15 after Rotation of " << phi * 180 / M_PI
                  << "degrees (+: counter-clockwise / -: clockwise)\n";
        for (uint k = 0; k < 15; k++) {
            X_temp = originXY[0] + (MP_X[k] - originXY[0]) * cos(phi) - (MP_Y[k] - originXY[1]) * sin(phi);
            Y_temp = originXY[1] + (MP_X[k] - originXY[0]) * sin(phi) + (MP_Y[k] - originXY[1]) * cos(phi);
            std::cout << "Name:  " << MP_name[k] << " ";
            std::cout << "\t\tX: " << X_temp << " ";
            std::cout << "\t\tY: " << Y_temp << " ";
            std::cout << "\t\tZ: " << MP_Z[k] << " " << std::endl;

            MP_X[k] = X_temp;
            MP_Y[k] = Y_temp;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Coordinates to Index Procedure
    // std::cout << "Converting Coordinates to Index..." << std::endl;
    for (uint k = 0; k < numberOfMeasurePoints; k++) {
        MP_k.push_back(
            gridBuilder->getGrid(MP_level[k])
                ->getSparseIndex(gridBuilder->getGrid(MP_level[k])->transCoordToIndex(MP_X[k], MP_Y[k], MP_Z[k])));
        if (MP_k[k] == -1) {
            std::cerr << "Error: Could not convert Coordinate to Sparse Index for MP " << k + 1 << std::endl;
        }
        std::cout << MP_name[k] << "\tID = "
                  << gridBuilder->getGrid(MP_level[k])
                         ->getSparseIndex(
                             gridBuilder->getGrid(MP_level[k])->transCoordToIndex(MP_X[k], MP_Y[k], MP_Z[k]))
                  << std::endl;
        // std::cout << "ID = " <<
        // gridBuilder->getGrid(0)->getSparseIndex(gridBuilder->getGrid(0)->transCoordToIndex(-0.500000,
        // -0.500000, 9.500000)) << std::endl;
    }
    // std::cout << "Done Converting Coordinates to Index..." << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Writing Procedure
    // std::cout << "Writing new file..." << std::endl;
    std::ofstream outFile(inputPath + "measurePoints.dat");

    outFile << numberOfMeasurePoints << std::endl;
    for (uint j = 0; j < numberOfMeasurePoints; j++) {
        outFile << MP_name[j] << " " << MP_k[j] << " " << MP_level[j] << std::endl;
        // std::cout << MP_name[j] << "\t" << MP_k[j] << "\t" << MP_level[j] << std::endl;
    }
    // std::cout << "Done writing new file..." << std::endl;
    outFile.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

std::string chooseVariation()
{
    switch (variant) {
        case 1:
            simulationName = "Var01_0deg_FD_a";
            break;
        case 2:
            simulationName = "Var02_0deg_FD_b";
            break;
        case 3:
            simulationName = "Var03_0deg_FD_c";
            break;
        case 4:
            simulationName = "Var04_0deg_SD_a";
            break;
        case 5:
            simulationName = "Var05_0deg_SD_b";
            break;
        case 6:
            simulationName = "Var06_0deg_SD_c";
            break;
        case 7:
            simulationName = "Var07_63deg_FD_a";
            break;
        case 8:
            simulationName = "Var08_63deg_FD_b";
            break;
        case 9:
            simulationName = "Var09_63deg_FD_c";
            break;
        case 10:
            simulationName = "Var10_63deg_SD_a";
            break;
        case 11:
            simulationName = "Var11_63deg_SD_b";
            break;
        case 12:
            simulationName = "Var12_63deg_SD_c";
            break;
        default:
            std::cerr << "Warning: no variant selected. Running with Default variant V01!" << std::endl;
            simulationName = "Var01_0deg_FD_a";
            rotationOfCity = 0.0;
            return simulationName;
    }

    if ((0 < variant) && (variant <= 6))
        rotationOfCity = (real)0.0;
    else if ((6 < variant) && (variant <= 12))
        rotationOfCity = (real)63.295;

    std::cout << "Variant selected. Simulation name is: " << simulationName << std::endl;
    std::cout << "Rotation selected: " << rotationOfCity << std::endl;

    return simulationName;
}


int main( int argc, char* argv[])
{
    try
    {
        vf::logging::Logger::initializeLogger();

        // assuming that the config files is stored parallel to this file.
        std::filesystem::path filePath = __FILE__;
        filePath.replace_filename("configDrivenCavity.txt");

        multipleLevel(filePath.string());
    }
    catch (const spdlog::spdlog_ex &ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    }
    catch (const std::bad_alloc& e)
    {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    }
    catch (const std::exception& e)
    {
        VF_LOG_CRITICAL("exception: {}", e.what());
    }
    catch (...)
    {
        VF_LOG_CRITICAL("Unknown exception!");
    }

   return 0;
}
