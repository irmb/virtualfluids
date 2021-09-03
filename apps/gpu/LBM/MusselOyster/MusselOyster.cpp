
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "basics/Core/DataTypes.h"
#include "basics/PointerDefinitions.h"
#include "basics/Core/VectorTypes.h"

#include "basics/Core/LbmOrGks.h"
#include "basics/Core/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"
#include "basics/Core/Logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridFactory.h"

#include "geometries/Sphere/Sphere.h"
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

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////

#include "utilities/communication.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string path("E:/temp/MusselOysterResults");
std::string gridPathParent = "E:/temp/GridMussel/";
std::string simulationName("MusselOyster");

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
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
	vf::gpu::Communicator* comm = vf::gpu::Communicator::getInstanz();
    vf::basics::ConfigurationFile config;
    std::cout << configPath << std::endl;
    config.load(configPath);
    SPtr<Parameter> para = std::make_shared<Parameter>(config, comm->getNummberOfProcess(), comm->getPID());



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = true;
    bool useMultiGPU = true;
    bool useStreams  = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::string bivalveType = "MUSSEL"; // "MUSSEL" "OYSTER"
    std::string gridPath(gridPathParent + bivalveType); // only for GridGenerator, for GridReader the gridPath needs to be set in the config file

    real dxGrid = (real)1.0;
    real vxLB = (real)0.051; // LB units
    real Re = (real)300.0;
    real viscosityLB = (vxLB * dxGrid) / Re;

    para->setVelocity(vxLB);
    para->setViscosity(viscosityLB);
    para->setVelocityRatio((real)58.82352941);
    para->setViscosityRatio((real)0.058823529);

    *logging::out << logging::Logger::INFO_HIGH << "velocity LB [dx/dt] = " << vxLB << " \n";
    *logging::out << logging::Logger::INFO_HIGH << "viscosity LB [dx^2/dt] = " << viscosityLB << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "velocity real [m/s] = " << vxLB * para->getVelocityRatio()<< " \n";
    *logging::out << logging::Logger::INFO_HIGH << "viscosity real [m^2/s] = " << viscosityLB * para->getViscosityRatio() << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "dxGrid = " << dxGrid << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "useGridGenerator = " << useGridGenerator << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "useMultiGPU = " << useMultiGPU << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "useStreams = " << useStreams << "\n";

    
    para->setTOut(10000);
    para->setTEnd(10000);

    para->setCalcDragLift(false);
    para->setUseWale(false);

    para->setOutputPath(path);
    para->setOutputPrefix(simulationName);
    para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());
    para->setPrintFiles(true);

    para->setMaxLevel(1);

    //para->setMainKernel("CumulantK17CompChim");
    if (useStreams)
        para->setUseStreams();
    para->setMainKernel("CumulantK17CompChimStream");
    *logging::out << logging::Logger::INFO_HIGH << "Kernel: " << para->getMainKernel() << "\n";

    if (useMultiGPU) {
        para->setDevices(std::vector<uint>{ (uint)0, (uint)1 });
        para->setMaxDev(2);
    } else 
        para->setDevices(std::vector<uint>{ (uint)0 });



    //////////////////////////////////////////////////////////////////////////


    if (useGridGenerator) {
        // bounding box mussel:
        const real bbxm = -18.0;
        const real bbxp = 58.0;
        const real bbym = -17.0;
        const real bbyp = 18.0;
        const real bbzm = -5.0;
        const real bbzp = 13.0;
        // bounding box oyster:
        // const real bbxm = 0.0;
        // const real bbxp = 115.0;
        // const real bbym = 0.0;
        // const real bbyp = 27.0;
        // const real bbzm = 0.0;
        // const real bbzp = 63.0;

        const real xGridMin  = bbxm - 40.0;
        const real xGridMax  = bbxp + 250.0;
        const real yGridMin  = bbym + 1.0;
        const real yGridMax  = bbyp + 60.0;
        const real zGridMin  = bbzm - 30.0;
        const real zGridMax  = bbzp + 30.0;

        TriangularMesh *bivalveSTL =
            TriangularMesh::make("C:/Users/Master/Documents/MasterAnna/STL/" + bivalveType + ".stl");
         //TriangularMesh* bivalveRef_1_STL =
         //    TriangularMesh::make("C:/Users/Master/Documents/MasterAnna/STL/" + bivalveType + "_Level1.stl");

        if (useMultiGPU) {
            const uint generatePart = vf::gpu::Communicator::getInstanz()->getPID();
            
            real overlap      = (real)8.0 * dxGrid;            
            const real ySplit = bbyp - 10.0;

            if (generatePart == 0) {
                gridBuilder->addCoarseGrid( xGridMin,   yGridMin,           zGridMin, 
                                            xGridMax,   ySplit+overlap,     zGridMax,   dxGrid);
            }
            if (generatePart == 1) {
                gridBuilder->addCoarseGrid(xGridMin,    ySplit-overlap,     zGridMin, 
                                           xGridMax,    yGridMax,           zGridMax,   dxGrid);
            }

             //gridBuilder->setNumberOfLayers(6, 8);
             //gridBuilder->addGrid(bivalveRef_1_STL, 1);

            gridBuilder->addGeometry(bivalveSTL);

            if (generatePart == 0)
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xGridMin,    xGridMax,
                                                                           yGridMin,    ySplit, 
                                                                           zGridMin,    zGridMax));
            if (generatePart == 1)
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xGridMin,    xGridMax, 
                                                                           ySplit,      yGridMax, 
                                                                           zGridMin,    zGridMax));
            
            gridBuilder->setPeriodicBoundaryCondition(false, false, true);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            if (generatePart == 0) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::PY, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PY, 1);
            }

            if (generatePart == 1) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::MY, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MY, 0);
            }

            //////////////////////////////////////////////////////////////////////////                       
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);
            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); 
            if (generatePart == 0)
                gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0); 
            if (generatePart == 1)
                gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);  
            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            //////////////////////////////////////////////////////////////////////////
            if (useStreams)
                gridBuilder->findFluidNodes(useStreams);

            //gridBuilder->writeGridsToVtk(path + "/" + bivalveType + "/grid/part" + std::to_string(generatePart) + "_");
            //gridBuilder->writeGridsToVtk(path + "/" + bivalveType + "/" + std::to_string(generatePart) + "/grid/");
            //gridBuilder->writeArrows(path + "/" + bivalveType + "/" + std::to_string(generatePart) + " /arrow");

            SimulationFileWriter::write(gridPath + "/" + std::to_string(generatePart) + "/", gridBuilder, FILEFORMAT::BINARY);
           
        } else {

            gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, xGridMax, yGridMax, zGridMax, dxGrid);

            //gridBuilder->setNumberOfLayers(6, 8);
            //gridBuilder->addGrid(bivalveRef_1_STL, 1);

            gridBuilder->addGeometry(bivalveSTL);

            gridBuilder->setPeriodicBoundaryCondition(false, false, true);

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////
            gridBuilder->setVelocityBoundaryCondition(SideType::PY, vxLB, 0.0, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
            gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
            gridBuilder->setVelocityBoundaryCondition(SideType::MX, vxLB, 0.0, 0.0);

            gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
            //////////////////////////////////////////////////////////////////////////
            if (useStreams)
                gridBuilder->findFluidNodes(useStreams);

            // gridBuilder->writeGridsToVtk("E:/temp/MusselOyster/" + bivalveType + "/grid/");
            // gridBuilder->writeArrows ("E:/temp/MusselOyster/" + bivalveType + "/arrow");

            SimulationFileWriter::write(gridPath, gridBuilder, FILEFORMAT::BINARY);
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // const real velocityLB = velocity * dt / dx; // LB units

        // const real vx = velocityLB / (real)sqrt(2.0); // LB units
        // const real vy = velocityLB / (real)sqrt(2.0); // LB units

        // const real viscosityLB = nx * velocityLB / Re; // LB units

        //*logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
        //*logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // para->setVelocity(velocityLB);
        // para->setViscosity(viscosityLB);

        // para->setVelocityRatio(velocity/ velocityLB);

        // para->setMainKernel("CumulantK17CompChim");

        // para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz)
        // {
        //          rho = (real)0.0;
        //          vx  = (real)0.0; //(6 * velocityLB * coordZ * (L - coordZ) / (L * L));
        //          vy  = (real)0.0;
        //          vz  = (real)0.0;
        //      });



       //return;
    }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

    SPtr<GridProvider> gridGenerator;
    if (useGridGenerator)
        gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);
    else {
        gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);
    }
           
    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    SPtr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
    SPtr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();
    sim.setFactories(kernelFactory, preProcessorFactory);
    sim.init(para, gridGenerator, fileWriter, cudaMemoryManager);

    //// Level 0
    //uint *geoSP   = para->getParH(0)->geoSP;
    //uint numGeoFluid = 0;
    //std::vector<int> sparseIndicesFluid;
    //for (uint i = 0; i < para->getParH(0)->size_Mat_SP; i++) {
    //    if (geoSP[i] == GEO_FLUID) {
    //        numGeoFluid++;
    //        sparseIndicesFluid.push_back(i);
    //    }
    //}
    //std::cout << ".....geoFluid level 0 " << numGeoFluid << ", num fluid nodes (new kernel)  " << para->getParH(0)->numberOfFluidNodes
    //          << std::endl;


    //for (uint i = 300000; i < 300003; i++)
    //    std::cout << ".....level 0: sparse index geoFluid \t" << sparseIndicesFluid[i] << ",    fluid nodes index  \t"
    //              << para->getParH(0)->fluidNodeIndices[i] << std::endl;

    //// Level 1
    //uint *geoSP1      = para->getParH(1)->geoSP;
    //uint numGeoFluid1 = 0;
    //std::vector<int> sparseIndicesFluid1;
    //for (uint i = 0; i < para->getParH(1)->size_Mat_SP; i++) {
    //    if (geoSP1[i] == GEO_FLUID) {
    //        numGeoFluid1++;
    //        sparseIndicesFluid1.push_back(i);
    //    }
    //}
    //std::cout << ".....geoFluid level 1 " << numGeoFluid1 << ", num fluid nodes (new kernel)  "
    //          << para->getParH(1)->numberOfFluidNodes << std::endl;

    //for (uint i = 300000; i < 300003; i++)
    //    std::cout << ".....level 1: sparse index geoFluid \t" << sparseIndicesFluid[i] << ",    fluid nodes index  \t"
    //              << para->getParH(0)->fluidNodeIndices[i] << std::endl;


    sim.run();
    sim.free();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
}

int main( int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    std::string str, str2; 
    if ( argv != NULL )
    {
        //str = static_cast<std::string>(argv[0]);
        
        try
        {
            //////////////////////////////////////////////////////////////////////////

			std::string targetPath;

			targetPath = __FILE__;

#ifdef _WIN32
            targetPath = targetPath.substr(0, targetPath.find_last_of('\\') + 1);
#else
            targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
#endif

			std::cout << targetPath << std::endl;

			multipleLevel(targetPath + "configMusselOyster.txt");

            //////////////////////////////////////////////////////////////////////////
		}
        catch (const std::bad_alloc& e)
        { 
            *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
        }
        catch (const std::exception& e)
        {   
            *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
        }
        catch (...)
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
        }
    }

   MPI_Finalize();
   return 0;
}
