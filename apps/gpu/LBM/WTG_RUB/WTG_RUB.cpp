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

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/LbmOrGks.h"
#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Core/Input/ConfigFileReader/ConfigFileReader.h"

#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridFactory.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"

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

const real L  = 1.0;

const real Re = 500.0;// 1000.0;

const real velocity  = 1.0;

int variant = 1;
real rotationOfCity;
real z_offset = 0.0; // only if baseplate is in use (currently not!! not important)

// 1: original setup of Lennard Lux (6 level, 4.0 cm -> 1.25 mm)
// 2: setup 1 of MSch               (4 level, 1.0 cm -> 1.25 mm)
// 3: setup 2 of MSch               (5 level, 1.6 cm -> 1.0  mm)
int setupDomain = 3;



std::string path("D:/out/WTG_RUB"); //Mollok
std::string inputPath("D:/out/WTG_RUB/input/");

std::string simulationName("RUB");

const uint timeStepStartOut = 0;
const uint timeStepOut = 10000;
const uint timeStepEnd = 100000;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void addFineGrids(SPtr<MultipleGridBuilder> gridBuilder, uint &maxLevel, real &rotationOfCity);

std::string chooseVariation();

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
	SPtr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance();

    std::cout << configPath << std::endl;
	SPtr<ConfigData> configData = configReader->readConfigFile(configPath.c_str());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real dx = 0;
    real viscosityLB = 1.0e-03;
    uint maxLevel    = 1;

    if (setupDomain == 1) {
        dx = 4;
        maxLevel    = 5;
        viscosityLB = 3.75e-06; // LB units
    } else if (setupDomain == 2) {
        dx = 1;   
        maxLevel    = 3;
        viscosityLB = 1.5e-05; // LB units
    } else if (setupDomain == 3) {
        dx          = 1.6;
        maxLevel    = 4;
        viscosityLB = 9.375e-06; // LB units
    }
    
    real x_min = 0.0;
    real x_max = 1250.0;
    real y_min = 0.0;
    real y_max = 190.0;
    real z_min = 0.0 + z_offset;
    real z_max = 160.0 + z_offset;

    //TriangularMesh *RubSTL      = TriangularMesh::make(inputPath + "stl/Var02_0deg_FD_b.stl");
    TriangularMesh *RubSTL      = TriangularMesh::make(inputPath + "stl/" + chooseVariation() + ".stl");
    // vector<real> originOfCityXY = { 600.0, y_max / 2, z_offset };


    gridBuilder->addCoarseGrid(x_min, y_min, z_min, 
                               x_max, y_max, z_max, dx);

    gridBuilder->setNumberOfLayers(0, 0);

    addFineGrids(gridBuilder, maxLevel, rotationOfCity);

    //// adding solid CityGeometry to gridbuilder
    gridBuilder->addGeometry(RubSTL);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<Parameter>    para         = Parameter::make(configData, comm);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const real velocityLB = 0.0844; // LB units

	//const real vx = velocityLB / (real)sqrt(2.0); // LB units
	//const real vy = velocityLB / (real)sqrt(2.0); // LB units

    *logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
    *logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	para->setDevices(std::vector<uint>{(uint)0});

    para->setOutputPath( path );
    para->setOutputPrefix( simulationName );

    para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());

    para->setPrintFiles(true);

    para->setMaxLevel(maxLevel);

    para->setVelocity(velocityLB);
    para->setViscosity(viscosityLB);

    para->setVelocityRatio(velocity/ velocityLB);

	para->setMainKernel("CumulantK17CompChim");

	para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = (real)0.0; //(6 * velocityLB * coordZ * (L - coordZ) / (L * L));
        vy  = (real)0.0;
        vz  = (real)0.0;
    });

    para->setTOut( timeStepOut );
    para->setTEnd( timeStepEnd );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);

    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    SPtr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
    SPtr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();
    sim.setFactories(kernelFactory, preProcessorFactory);
    sim.init(para, gridGenerator, fileWriter, cudaMemoryManager);
    sim.run();
    sim.free();

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
        auto FG1 = new Cuboid(-20, -20, -5 + z_offset, 800, 200, 75 + z_offset);

        // FineGrid Level 2 -> dx = 1 cm; lvl 2
        auto FG2_1 = new Cuboid(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG2_2 = new Cuboid(500, -20,  5 + z_offset, 680, 210, 50 + z_offset);
        auto FG2   = new Conglomerate();
        FG2->add(FG2_1);
        FG2->add(FG2_2);

        // FineGrid Level 3 ->dx = 5 mm; lvl 3
        auto FG3_1 = new Cuboid(517, -20, -5 + z_offset, 665, 200, 30 + z_offset);
        auto FG3_2 = new Cuboid(550, 58, -5 + z_offset, 650, 132, 40 + z_offset);
        auto FG3   = new Conglomerate();
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
                            TriangularMesh *FG4 = TriangularMesh::make(inputPath + "stl/FG4_0deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        } else {
                            TriangularMesh *FG4 = TriangularMesh::make(inputPath + "stl/FG4_63deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        }

                        if (maxLevel == 5) {
                            if (rotationOfCity == 0.0) {
                                TriangularMesh *FG5 = TriangularMesh::make(inputPath + "stl/FG5_0deg.stl");
                                gridBuilder->addGrid(FG5, 5);
                            } else {
                                TriangularMesh *FG5 = TriangularMesh::make(inputPath + "stl/FG5_63deg.stl");
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
        //auto FG1_1 = new Cuboid(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG1_1 = new Cuboid(-20, -20, -5 + z_offset, 760, 200, 20 + z_offset);
        auto FG1_2 = new Cuboid(500, -20,  5 + z_offset, 680, 210, 50 + z_offset);
        auto FG1   = new Conglomerate();
        FG1->add(FG1_1);
        FG1->add(FG1_2);

        // Adding FineGrids 1 to 5 depending on maxLevel, FG4 and FG5 require different STL-files depending on
        // rotationOfCity
        if (maxLevel >= 1) {
            gridBuilder->addGrid(FG1, 1);
            if (maxLevel >= 2) {
                if (rotationOfCity == 0.0) {
                    TriangularMesh *FG2 = TriangularMesh::make(inputPath + "stl/FG4_0deg.stl");
                    gridBuilder->addGrid(FG2, 2);
                } else {
                    TriangularMesh *FG2 = TriangularMesh::make(inputPath + "stl/FG4_63deg.stl");
                    gridBuilder->addGrid(FG2, 2);
                }

                if (maxLevel == 3) {
                    if (rotationOfCity == 0.0) {
                        TriangularMesh *FG3 = TriangularMesh::make(inputPath + "stl/FG5_0deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    } else {
                        TriangularMesh *FG3 = TriangularMesh::make(inputPath + "stl/FG5_63deg.stl");
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
        //auto FG1 = new Cuboid(-20, -20, -5 + z_offset, 800, 200, 75 + z_offset);

        // FineGrid Level 1 -> dx = 8.0 mm; lvl 1
        auto FG1_1 = new Cuboid(-20, -20, -5 + z_offset, 780, 200, 30 + z_offset);
        auto FG1_2 = new Cuboid(500, -20,  5 + z_offset, 720, 210, 75 + z_offset);
        auto FG1   = new Conglomerate();
        FG1->add(FG1_1);
        FG1->add(FG1_2);

        // FineGrid Level 2 -> dx = 4.0 mm; lvl 2
        auto FG2_1 = new Cuboid(-20, -20, -5 + z_offset, 760, 200, 10 + z_offset);
        auto FG2_2 = new Cuboid(520, -20,  5 + z_offset, 700, 210, 50 + z_offset);
        auto FG2   = new Conglomerate();
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
                        TriangularMesh *FG3 = TriangularMesh::make(inputPath + "stl/FG4_0deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    } else {
                        TriangularMesh *FG3 = TriangularMesh::make(inputPath + "stl/FG4_63deg.stl");
                        gridBuilder->addGrid(FG3, 3);
                    }

                    if (maxLevel == 4) {
                        if (rotationOfCity == 0.0) {
                            TriangularMesh *FG4 = TriangularMesh::make(inputPath + "stl/FG5_0deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        } else {
                            TriangularMesh *FG4 = TriangularMesh::make(inputPath + "stl/FG5_63deg.stl");
                            gridBuilder->addGrid(FG4, 4);
                        }
                    }
                }
            }
        }
    }



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
        rotationOfCity = 0.0;
    else if ((6 < variant) && (variant <= 12))
        rotationOfCity = 63.295;

    std::cout << "Variant selected. Simulation name is: " << simulationName << std::endl;
    std::cout << "Rotation selected: " << rotationOfCity << std::endl;

    return simulationName;
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

			targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);



			std::cout << targetPath << std::endl;

			multipleLevel(targetPath + "configDrivenCavity.txt");

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
