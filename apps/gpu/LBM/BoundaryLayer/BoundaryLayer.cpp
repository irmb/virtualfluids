
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>
#include <numeric>

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/StringUtilities/StringUtil.h"

#include "Core/VectorTypes.h"

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>


//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridFactory.h"

#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/ActuatorLine.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlanarAverageProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/WallModelProbe.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/TurbulenceModels/TurbulenceModelFactory.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

#include "utilities/communication.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string path(".");

std::string simulationName("BoundayLayer");

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
    
    const int  nProcs = communicator.getNummberOfProcess();
    const uint procID = vf::gpu::Communicator::getInstance().getPID();
    std::vector<uint> devices(10);
    std::iota(devices.begin(), devices.end(), 0);
    para->setDevices(devices);
    para->setMaxDev(nProcs);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //          U s e r    s e t t i n g s
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    LbmOrGks lbmOrGks = LBM;


    real H = 1000.0; // boundary layer height in m
    const real L_x = config.contains("L_x")? config.getValue<real>("L_x"): 6*H;
    const real L_y = config.contains("L_y")? config.getValue<real>("L_y"): 4*H;
    const real L_z = config.contains("L_z")? config.getValue<real>("L_z"): 1*H;
    H = L_z;

    const real z0  = 0.1; // roughness length in m
    const real u_star = 0.4; //friction velocity in m/s
    const real kappa = 0.4; // von Karman constant

    const real viscosity = 1.56e-5;

    const real velocity  = 0.5*u_star/kappa*log(L_z/z0); //0.5 times max mean velocity at the top in m/s

    const real mach = config.contains("Ma")? config.getValue<real>("Ma"): 0.1;

    const uint nodes_per_H = config.contains("nz")? config.getValue<uint>("nz"): 64;

    // all in s
    const float tStartOut   = config.getValue<real>("tStartOut");
    const float tOut        = config.getValue<real>("tOut");
    const float tEnd        = config.getValue<real>("tEnd"); // total time of simulation

    const float tStartAveraging     =  config.getValue<real>("tStartAveraging");
    const float tStartTmpAveraging  =  config.getValue<real>("tStartTmpAveraging");
    const float tAveraging          =  config.getValue<real>("tAveraging");
    const float tStartOutProbe      =  config.getValue<real>("tStartOutProbe");
    const float tOutProbe           =  config.getValue<real>("tOutProbe");


    const real dx = L_z/real(nodes_per_H);

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

    bool useStreams = (nProcs > 1 ? true: false);
    // useStreams=true;
    para->setUseStreams(useStreams);
    para->setMainKernel("CumulantK17Almighty");
    para->setIsBodyForce( config.getValue<bool>("bodyForce") );

    para->setTimestepStartOut(uint(tStartOut/dt) );
    para->setTimestepOut( uint(tOut/dt) );
    para->setTimestepEnd( uint(tEnd/dt) );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<TurbulenceModelFactory> tmFactory = SPtr<TurbulenceModelFactory>( new TurbulenceModelFactory(para) );
    tmFactory->readConfigFile( config );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real xSplit = L_x/nProcs;
    const real overlap = 8.0*dx;

    real xMin      =  procID    * xSplit;
    real xMax      = (procID+1) * xSplit;
    real xGridMin  =  procID    * xSplit;
    real xGridMax  = (procID+1) * xSplit;

    real yMin      = 0.0;
    real yMax      = L_y;
    real zMin      = 0.0;
    real zMax      = L_z; 

    bool isFirstSubDomain = (procID == 0        && nProcs > 1)?                    true: false;
    bool isLastSubDomain  = (procID == nProcs-1 && nProcs > 1)?                    true: false;
    bool isMidSubDomain   = (!isFirstSubDomain && !isLastSubDomain && nProcs > 1)? true: false;
    
    if(isFirstSubDomain || isMidSubDomain)
    {
        xGridMax += overlap;
        xGridMin -= overlap;
    }
    if(isLastSubDomain || isMidSubDomain)
    {
        xGridMax += overlap;
        xGridMin -= overlap;
    }

    gridBuilder->addCoarseGrid( xGridMin,  0.0,  0.0,
                                xGridMax,  L_y,  L_z, dx);
    if(false)
    {
        gridBuilder->setNumberOfLayers(12, 8);
        gridBuilder->addGrid( new Cuboid( 0.0, 0.0, 0.0, L_x,  L_y,  0.3*L_z) , 1 );
        para->setMaxLevel(2);
    }

    if(nProcs > 1)
    {
            gridBuilder->setSubDomainBox(
                        std::make_shared<BoundingBox>(xMin, xMax, yMin, yMax, zMin, zMax));        
            gridBuilder->setPeriodicBoundaryCondition(false, true, false);
    }
    else         
    { 
        gridBuilder->setPeriodicBoundaryCondition(true, true, false);
    }

	gridBuilder->buildGrids(lbmOrGks, true); // buildGrids() has to be called before setting the BCs!!!!

    std::cout << "nProcs: "<< nProcs << "Proc: " << procID << " isFirstSubDomain: " << isFirstSubDomain << " isLastSubDomain: " << isLastSubDomain << " isMidSubDomain: " << isMidSubDomain << std::endl;
    
    if(nProcs > 1){
        if (isFirstSubDomain || isMidSubDomain) {
            gridBuilder->findCommunicationIndices(CommunicationDirections::PX, lbmOrGks);
            gridBuilder->setCommunicationProcess(CommunicationDirections::PX, procID+1);
        }

        if (isLastSubDomain || isMidSubDomain) {
            gridBuilder->findCommunicationIndices(CommunicationDirections::MX, lbmOrGks);
            gridBuilder->setCommunicationProcess(CommunicationDirections::MX, procID-1);
        }

        if (isFirstSubDomain) {
            gridBuilder->findCommunicationIndices(CommunicationDirections::MX, lbmOrGks);
            gridBuilder->setCommunicationProcess(CommunicationDirections::MX, nProcs-1);
        }

        if (isLastSubDomain) {
            gridBuilder->findCommunicationIndices(CommunicationDirections::PX, lbmOrGks);
            gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 0);
        }
    }
    uint samplingOffset = 2;
    gridBuilder->setStressBoundaryCondition(SideType::MZ,
                                            0.0, 0.0, 1.0,              // wall normals
                                            samplingOffset, z0/dx);     // wall model settinng
    para->setHasWallModelMonitor(true);
    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressPressureBounceBack);

    gridBuilder->setSlipBoundaryCondition(SideType::PZ,  0.0,  0.0, 0.0);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipBounceBack); 
    

    real cPi = 3.1415926535897932384626433832795;
    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = (u_star/0.4 * log(coordZ/z0) + 2.0*sin(cPi*16.0f*coordX/L_x)*sin(cPi*8.0f*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1))  * dt / dx; 
        vy  = 2.0*sin(cPi*16.0f*coordX/L_x)*sin(cPi*8.0f*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1)  * dt / dx; 
        vz  = 8.0*u_star/0.4*(sin(cPi*8.0*coordY/H)*sin(cPi*8.0*coordZ/H)+sin(cPi*8.0*coordX/L_x))/(pow(L_z/2.0-coordZ, c2o1)+c1o1) * dt / dx;
    });
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(isFirstSubDomain || nProcs == 1)
    {
        SPtr<PlanarAverageProbe> planarAverageProbe = SPtr<PlanarAverageProbe>( new PlanarAverageProbe("planeProbe", para->getOutputPath(), tStartAveraging/dt, tStartTmpAveraging/dt, tAveraging/dt , tStartOutProbe/dt, tOutProbe/dt, 'z') );
        planarAverageProbe->addAllAvailableStatistics();
        planarAverageProbe->setFileNameToNOut();
        para->addProbe( planarAverageProbe );

        para->setHasWallModelMonitor(true);
        SPtr<WallModelProbe> wallModelProbe = SPtr<WallModelProbe>( new WallModelProbe("wallModelProbe", para->getOutputPath(), tStartAveraging/dt, tStartTmpAveraging/dt, tAveraging/dt/4.0 , tStartOutProbe/dt, tOutProbe/dt) );
        wallModelProbe->addAllAvailableStatistics();
        wallModelProbe->setFileNameToNOut();
        wallModelProbe->setForceOutputToStress(true);
        if(para->getIsBodyForce())
            wallModelProbe->setEvaluatePressureGradient(true);
        para->addProbe( wallModelProbe );

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
