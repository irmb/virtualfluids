
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/StringUtilities/StringUtil.h"

#include "Core/VectorTypes.h"

#include <basics/config/ConfigurationFile.h>
#include "lbm/constants/NumericConstants.h"

#include <logger/Logger.h>


//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"

#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/VelocitySetter/VelocitySetter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlanarAverageProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/WallModelProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/PrecursorWriter.h"
#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/TurbulenceModels/TurbulenceModelFactory.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string path(".");

std::string simulationName("BoundaryLayer");

using namespace vf::lbm::constant;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void multipleLevel(const std::string& configPath)
{

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    auto gridFactory = GridFactory::make();
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();

    vf::basics::ConfigurationFile config;
    config.load(configPath);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNummberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //          U s e r    s e t t i n g s
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    LbmOrGks lbmOrGks = LBM;

    const real H = config.getValue("boundaryLayerHeight", 1000.0); // boundary layer height in m

    const real L_x = 6*H;
    const real L_y = 4*H;
    const real L_z = H;

    const real z0  = config.getValue("z0", 0.1f); // roughness length in m
    const real u_star = config.getValue("u_star", 0.4f); //friction velocity in m/s
    const real kappa = config.getValue("vonKarmanConstant", 0.4f); // von Karman constant

    const real viscosity = config.getValue("viscosity", 1.56e-5f);

    const real velocity  = 0.5f*u_star/kappa*log(H/z0+1.f); //0.5 times max mean velocity at the top in m/s

    const real mach = config.getValue<real>("Ma", 0.1);

    const uint nodes_per_H = config.getValue<uint>("nz", 64);

    const bool writePrecursor = config.getValue("writePrecursor", false);
    int nTWritePrecursor; real tStartPrecursor, posXPrecursor;
    if(writePrecursor)
    {
        nTWritePrecursor      = config.getValue<int>("nTimestepsWritePrecursor");
        tStartPrecursor      = config.getValue<real>("tStartPrecursor");
        posXPrecursor        = config.getValue<real>("posXPrecursor");
    }

    const bool readPrecursor = config.getValue("readPrecursor", false);
    int nTReadPrecursor;
    std::string precursorFile;
    if(readPrecursor)
    {
        nTReadPrecursor = config.getValue<int>("nTimestepsReadPrecursor");
        precursorFile = config.getValue<std::string>("precursorFile");
    }
    // all in s
    const float tStartOut   = config.getValue<real>("tStartOut");
    const float tOut        = config.getValue<real>("tOut");
    const float tEnd        = config.getValue<real>("tEnd"); // total time of simulation

    const float tStartAveraging     =  config.getValue<real>("tStartAveraging");
    const float tStartTmpAveraging  =  config.getValue<real>("tStartTmpAveraging");
    const float tAveraging          =  config.getValue<real>("tAveraging");
    const float tStartOutProbe      =  config.getValue<real>("tStartOutProbe");
    const float tOutProbe           =  config.getValue<real>("tOutProbe");


    const real dx = H/real(nodes_per_H);

    const real dt = dx * mach / (sqrt(3) * velocity);

    const real velocityLB = velocity * dt / dx; // LB units

    const real viscosityLB = viscosity * dt / (dx * dx); // LB units

    const real pressureGradient = u_star * u_star / H ;
    const real pressureGradientLB = pressureGradient * (dt*dt)/dx; // LB units

    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("dt   = {}", dt);
    VF_LOG_INFO("dx   = {}", dx);
    VF_LOG_INFO("viscosity [10^8 dx^2/dt] = {}", viscosityLB*1e8);
    VF_LOG_INFO("u* /(dx/dt) = {}", u_star*dt/dx);
    VF_LOG_INFO("dpdx  = {}", pressureGradient);
    VF_LOG_INFO("dpdx /(dx/dt^2) = {}", pressureGradientLB);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setOutputPrefix( simulationName );

    para->setPrintFiles(true);

    para->setForcing(pressureGradientLB, 0, 0);
    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio( dx / dt );
    para->setViscosityRatio( dx*dx/dt );
    para->setDensityRatio( 1.0 );

    para->setMainKernel("TurbulentViscosityCumulantK17CompChim");

    para->setIsBodyForce( config.getValue<bool>("bodyForce") );

    para->setTimestepStartOut(uint(tStartOut/dt) );
    para->setTimestepOut( uint(tOut/dt) );
    para->setTimestepEnd( uint(tEnd/dt) );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<TurbulenceModelFactory> tmFactory = SPtr<TurbulenceModelFactory>( new TurbulenceModelFactory(para) );
    tmFactory->readConfigFile( config );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0,
                                L_x,  L_y,  L_z, dx);
    // gridBuilder->setNumberOfLayers(12, 8);

    // gridBuilder->addGrid( new Cuboid( 0.0, 0.0, 0.0, L_x,  L_y,  0.3*L_z) , 1 );
    // para->setMaxLevel(2);

    gridBuilder->setPeriodicBoundaryCondition(!readPrecursor, true, false);

	gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

    
    if(readPrecursor)
    {

        uint samplingOffset = 2;
        gridBuilder->setStressBoundaryCondition(SideType::MZ,
                                            0.0, 0.0, 1.0,              // wall normals
                                            samplingOffset, z0/dx);     // wall model settinng
        para->setHasWallModelMonitor(true);
        auto precursor = createFileCollection(precursorFile, FileType::VTK);
        
        gridBuilder->setSlipBoundaryCondition(SideType::PZ,  0.0f,  0.0f, -1.0f);

        gridBuilder->setPrecursorBoundaryCondition(SideType::MX, precursor, nTReadPrecursor);
        
        gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.f);
    } 
    else
    {
        gridBuilder->setSlipBoundaryCondition(SideType::PZ,  0.0,  0.0, -1.0);

        uint samplingOffset = 2;
        gridBuilder->setStressBoundaryCondition(SideType::MZ,
                                            0.0, 0.0, 1.0,              // wall normals
                                            samplingOffset, z0/dx);     // wall model settinng
        para->setHasWallModelMonitor(true);
    }



    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressPressureBounceBack);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipBounceBack); 
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflectivePressureCorrection);
    para->setOutflowPressureCorrectionFactor(0.0); 
    




    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = (u_star/0.4 * log(coordZ/z0) + 2.0*sin(cPi*16.0f*coordX/L_x)*sin(cPi*8.0f*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1))  * dt / dx; 
        vy  = 2.0*sin(cPi*16.0f*coordX/L_x)*sin(cPi*8.0f*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1)  * dt / dx; 
        vz  = 8.0*u_star/0.4*(sin(cPi*8.0*coordY/H)*sin(cPi*8.0*coordZ/H)+sin(cPi*8.0*coordX/L_x))/(pow(L_z/2.0-coordZ, c2o1)+c1o1) * dt / dx;
        if(readPrecursor)
        {
            vx  = u_star/0.4 * log(coordZ/z0)  * dt / dx; 
            vy  = c0o1; 
            vz  = c0o1;
        }
    });


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<PlanarAverageProbe> planarAverageProbe = SPtr<PlanarAverageProbe>( new PlanarAverageProbe("horizontalPlanes", para->getOutputPath(), 0, tStartTmpAveraging/dt, tAveraging/dt , tStartOutProbe/dt, tOutProbe/dt, 'z') );
    planarAverageProbe->addAllAvailableStatistics();
    planarAverageProbe->setFileNameToNOut();
    para->addProbe( planarAverageProbe );

    SPtr<WallModelProbe> wallModelProbe = SPtr<WallModelProbe>( new WallModelProbe("wallModelProbe", para->getOutputPath(), 0, tStartTmpAveraging/dt, tAveraging/dt/4.0 , tStartOutProbe/dt, tOutProbe/dt) );
    wallModelProbe->addAllAvailableStatistics();
    wallModelProbe->setFileNameToNOut();
    wallModelProbe->setForceOutputToStress(true);
    if(para->getIsBodyForce())
        wallModelProbe->setEvaluatePressureGradient(true);
    para->addProbe( wallModelProbe );

    SPtr<PlaneProbe> planeProbe1 = SPtr<PlaneProbe>( new PlaneProbe("planeProbe_1", para->getOutputPath(), tStartAveraging/dt, 10, tStartOutProbe/dt, tOutProbe/dt) );
    planeProbe1->setProbePlane(100.0, 0.0, 0, dx, L_y, L_z);
    planeProbe1->addAllAvailableStatistics();
    para->addProbe( planeProbe1 );

    if(readPrecursor)
    {
        SPtr<PlaneProbe> planeProbe2 = SPtr<PlaneProbe>( new PlaneProbe("planeProbe_2", para->getOutputPath(), tStartAveraging/dt, 10, tStartOutProbe/dt, tOutProbe/dt) );
        planeProbe2->setProbePlane(1000.0, 0.0, 0, dx, L_y, L_z);
        planeProbe2->addAllAvailableStatistics();
        para->addProbe( planeProbe2 );

        SPtr<PlaneProbe> planeProbe3 = SPtr<PlaneProbe>( new PlaneProbe("planeProbe_3", para->getOutputPath(), tStartAveraging/dt, 10, tStartOutProbe/dt, tOutProbe/dt) );
        planeProbe3->setProbePlane(1500.0, 0.0, 0, dx, L_y, L_z);
        planeProbe3->addAllAvailableStatistics();
        para->addProbe( planeProbe3 );

        SPtr<PlaneProbe> planeProbe4 = SPtr<PlaneProbe>( new PlaneProbe("planeProbe_4", para->getOutputPath(), tStartAveraging/dt, 10, tStartOutProbe/dt, tOutProbe/dt) );
        planeProbe4->setProbePlane(2000.0, 0.0, 0, dx, L_y, L_z);
        planeProbe4->addAllAvailableStatistics();
        para->addProbe( planeProbe4 );

        SPtr<PlaneProbe> planeProbe5 = SPtr<PlaneProbe>( new PlaneProbe("planeProbe_5", para->getOutputPath(), tStartAveraging/dt, 10, tStartOutProbe/dt, tOutProbe/dt) );
        planeProbe5->setProbePlane(2500.0, 0.0, 0, dx, L_y, L_z);
        planeProbe5->addAllAvailableStatistics();
        para->addProbe( planeProbe5 );

        SPtr<PlaneProbe> planeProbe6 = SPtr<PlaneProbe>( new PlaneProbe("planeProbe_6", para->getOutputPath(), tStartAveraging/dt, 10, tStartOutProbe/dt, tOutProbe/dt) );
        planeProbe6->setProbePlane(0.0, L_y/2.0, 0, L_x, dx, L_z);
        planeProbe6->addAllAvailableStatistics();
        para->addProbe( planeProbe6 );
    }


    if(writePrecursor)
    {
        SPtr<PrecursorWriter> precursorWriter = SPtr<PrecursorWriter>( new PrecursorWriter("precursorWriter", para->getOutputPath()+"/precursor", posXPrecursor, 0, L_y, 0, L_z, tStartPrecursor/dt, nTWritePrecursor, OutputVariable::Distributions) );
        para->addProbe(precursorWriter);
    }

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
    auto gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, tmFactory);
    sim.run();
}

int main( int argc, char* argv[])
{
    if ( argv != NULL )
    {
        try
        {
            vf::logging::Logger::initalizeLogger();

            if( argc > 1){ path = argv[1]; }

            multipleLevel(path + "/configBoundaryLayer.txt");
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
    }
    return 0;
}
