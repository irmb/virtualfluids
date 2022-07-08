
#define _USE_MATH_DEFINES
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/StringUtilities/StringUtil.h"

#include "Core/VectorTypes.h"

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

#include "GridGenerator/geometries/Sphere/Sphere.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const real L = 1.0;
const real dSphere = 0.2;
const real Re = 1000.0; // related to the sphere's diameter
const real velocity = 1.0;
const real dt = (real)0.5e-3;
const uint nx = 64;

const uint timeStepOut = 10000;
const uint timeStepEnd = 100000;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string &configPath)
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    vf::gpu::Communicator &communicator = vf::gpu::Communicator::getInstance();

    auto gridFactory = GridFactory::make();
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dx = L / real(nx);

    gridBuilder->addCoarseGrid(-1.0 * L, -0.5 * L, -0.5 * L, 
                                6.0 * L,  0.5 * L,  0.5 * L, dx);

    Object *sphere = new Sphere(0.0, 0.0, 0.0, dSphere / 2.0);
    gridBuilder->addGeometry(sphere);

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(LBM, false); // buildGrids() has to be called before setting the BCs!!!!

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vf::basics::ConfigurationFile config;
    config.load(configPath);

    SPtr<Parameter> para =
        std::make_shared<Parameter>(config, communicator.getNummberOfProcess(), communicator.getPID());
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real velocityLB = velocity * dt / dx; // LB units

    const real viscosityLB = (dSphere / dx) * velocityLB / Re; // LB units

    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("viscosity [dx^2/dt] = {}", viscosityLB);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setDevices(std::vector<uint>{ (uint)0 });

    para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());

    para->setPrintFiles(true);

    para->setMaxLevel(2);

    para->setVelocity(velocityLB);
    para->setViscosity(viscosityLB);

    para->setVelocityRatio(velocity / velocityLB);
    para->setDensityRatio((real)1.0);

    // para->setMainKernel("CumulantK17CompChim");

    // para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
    //     rho = (real)0.0;
    //     vx =  (real)velocityLB;
    //     vy =  (real)0.0;
    //     vz =  (real)0.0;
    // });

    para->setTOut(timeStepOut);
    para->setTEnd(timeStepEnd);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);

    gridBuilder->setSlipBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::PZ, 0.0, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);

    // gridBuilder->setNoSlipBoundaryCondition(SideType::GEOMETRY); // not working yet, use veloBC
    gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, velocityLB * 2.0, 0.0);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure boundary condition last

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
    // bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // gridBuilder->writeGridsToVtk("grid/");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    auto gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
    sim.run();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

int main(int argc, char *argv[])
{
    try {
        vf::logging::Logger::initalizeLogger();

        // assuming that the config files is stored parallel to this file.
        std::filesystem::path filePath = __FILE__;
        filePath.replace_filename("config.txt");

        multipleLevel(filePath.string());
    } catch (const spdlog::spdlog_ex &ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    } catch (const std::bad_alloc &e) {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    } catch (const std::exception &e) {
        VF_LOG_CRITICAL("exception: {}", e.what());
    } catch (...) {
        VF_LOG_CRITICAL("Unknown exception!");
    }

    return 0;
}
