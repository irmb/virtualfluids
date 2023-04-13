//#define MPI_LOGGING


#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <iostream>
#include <stdexcept>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "metis.h"

#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "global.h"

#include "geometries/Sphere/Sphere.h"
#include "geometries/VerticalCylinder/VerticalCylinder.h"
#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMesh.h"
#include "geometries/Conglomerate/Conglomerate.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/GridFactory.h"

#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"

#include "utilities/math/Math.h"
#include "utilities/communication.h"
#include "utilities/transformator/TransformatorImp.h"

using namespace vf::gpu;

std::string getGridPath(std::shared_ptr<Parameter> para, std::string Gridpath)
{
    if (para->getNumprocs() == 1)
        return Gridpath + "/";
    
    return Gridpath + "/" + StringUtil::toString(para->getMyID()) + "/";
}

void setParameters(std::shared_ptr<Parameter> para, std::unique_ptr<input::Input> &input)
{
	Communicator* comm = Communicator::getInstanz();

	para->setMaxDev(StringUtil::toInt(input->getValue("NumberOfDevices")));
	para->setNumprocs(comm->getNumberOfProcess());
	para->setDevices(StringUtil::toUintVector(input->getValue("Devices")));
	para->setMyID(comm->getPID());
	
	std::string _path = input->getValue("Path");
    std::string _prefix = input->getValue("Prefix");
    std::string _gridpath = input->getValue("GridPath");
    std::string gridPath = getGridPath(para, _gridpath);
    para->setOutputPath(_path);
    para->setOutputPrefix(_prefix);
    para->setPathAndFilename(_path + "/" + _prefix);
    para->setPrintFiles(false);
    para->setPrintFiles(StringUtil::toBool(input->getValue("WriteGrid")));
    para->setGeometryValues(StringUtil::toBool(input->getValue("GeometryValues")));
    para->setCalc2ndOrderMoments(StringUtil::toBool(input->getValue("calc2ndOrderMoments")));
    para->setCalc3rdOrderMoments(StringUtil::toBool(input->getValue("calc3rdOrderMoments")));
    para->setCalcHighOrderMoments(StringUtil::toBool(input->getValue("calcHigherOrderMoments")));
    para->setReadGeo(StringUtil::toBool(input->getValue("ReadGeometry")));
    para->setCalcMedian(StringUtil::toBool(input->getValue("calcMedian")));
    para->setConcFile(StringUtil::toBool(input->getValue("UseConcFile")));
    para->setUseMeasurePoints(StringUtil::toBool(input->getValue("UseMeasurePoints")));
    para->setUseWale(StringUtil::toBool(input->getValue("UseWale")));
    para->setSimulatePorousMedia(StringUtil::toBool(input->getValue("SimulatePorousMedia")));
    para->setD3Qxx(StringUtil::toInt(input->getValue("D3Qxx")));
    para->setTimestepEnd(StringUtil::toInt(input->getValue("TimeEnd")));
    para->setTimestepOut(StringUtil::toInt(input->getValue("TimeOut")));
    para->setTimestepStartOut(StringUtil::toInt(input->getValue("TimeStartOut")));
    para->setTimeCalcMedStart(StringUtil::toInt(input->getValue("TimeStartCalcMedian")));
    para->setTimeCalcMedEnd(StringUtil::toInt(input->getValue("TimeEndCalcMedian")));
    para->setPressInID(StringUtil::toInt(input->getValue("PressInID")));
    para->setPressOutID(StringUtil::toInt(input->getValue("PressOutID")));
    para->setPressInZ(StringUtil::toInt(input->getValue("PressInZ")));
    para->setPressOutZ(StringUtil::toInt(input->getValue("PressOutZ")));
    //////////////////////////////////////////////////////////////////////////
    para->setDiffOn(StringUtil::toBool(input->getValue("DiffOn")));
    para->setDiffMod(StringUtil::toInt(input->getValue("DiffMod")));
    para->setDiffusivity(StringUtil::toFloat(input->getValue("Diffusivity")));
    para->setTemperatureInit(StringUtil::toFloat(input->getValue("Temp")));
    para->setTemperatureBC(StringUtil::toFloat(input->getValue("TempBC")));
    //////////////////////////////////////////////////////////////////////////
    para->setViscosityLB(StringUtil::toFloat(input->getValue("Viscosity_LB")));
    para->setVelocityLB(StringUtil::toFloat(input->getValue("Velocity_LB")));
    para->setViscosityRatio(StringUtil::toFloat(input->getValue("Viscosity_Ratio_World_to_LB")));
    para->setVelocityRatio(StringUtil::toFloat(input->getValue("Velocity_Ratio_World_to_LB")));
    para->setDensityRatio(StringUtil::toFloat(input->getValue("Density_Ratio_World_to_LB")));
    para->setPressRatio(StringUtil::toFloat(input->getValue("Delta_Press")));
    para->setRealX(StringUtil::toFloat(input->getValue("SliceRealX")));
    para->setRealY(StringUtil::toFloat(input->getValue("SliceRealY")));
    para->setFactorPressBC(StringUtil::toFloat(input->getValue("dfpbc")));
    para->setGeometryFileC(input->getValue("GeometryC"));
    para->setGeometryFileM(input->getValue("GeometryM"));
    para->setGeometryFileF(input->getValue("GeometryF"));
    //////////////////////////////////////////////////////////////////////////
    para->setgeoVec(gridPath + input->getValue("geoVec"));
    para->setcoordX(gridPath + input->getValue("coordX"));
    para->setcoordY(gridPath + input->getValue("coordY"));
    para->setcoordZ(gridPath + input->getValue("coordZ"));
    para->setneighborX(gridPath + input->getValue("neighborX"));
    para->setneighborY(gridPath + input->getValue("neighborY"));
    para->setneighborZ(gridPath + input->getValue("neighborZ"));
    para->setscaleCFC(gridPath + input->getValue("scaleCFC"));
    para->setscaleCFF(gridPath + input->getValue("scaleCFF"));
    para->setscaleFCC(gridPath + input->getValue("scaleFCC"));
    para->setscaleFCF(gridPath + input->getValue("scaleFCF"));
    para->setscaleOffsetCF(gridPath + input->getValue("scaleOffsetCF"));
    para->setscaleOffsetFC(gridPath + input->getValue("scaleOffsetFC"));
    para->setgeomBoundaryBcQs(gridPath + input->getValue("geomBoundaryBcQs"));
    para->setgeomBoundaryBcValues(gridPath + input->getValue("geomBoundaryBcValues"));
    para->setinletBcQs(gridPath + input->getValue("inletBcQs"));
    para->setinletBcValues(gridPath + input->getValue("inletBcValues"));
    para->setoutletBcQs(gridPath + input->getValue("outletBcQs"));
    para->setoutletBcValues(gridPath + input->getValue("outletBcValues"));
    para->settopBcQs(gridPath + input->getValue("topBcQs"));
    para->settopBcValues(gridPath + input->getValue("topBcValues"));
    para->setbottomBcQs(gridPath + input->getValue("bottomBcQs"));
    para->setbottomBcValues(gridPath + input->getValue("bottomBcValues"));
    para->setfrontBcQs(gridPath + input->getValue("frontBcQs"));
    para->setfrontBcValues(gridPath + input->getValue("frontBcValues"));
    para->setbackBcQs(gridPath + input->getValue("backBcQs"));
    para->setbackBcValues(gridPath + input->getValue("backBcValues"));
    para->setnumberNodes(gridPath + input->getValue("numberNodes"));
    para->setLBMvsSI(gridPath + input->getValue("LBMvsSI"));
    //////////////////////////////gridPath + ////////////////////////////////////////////
    para->setmeasurePoints(gridPath + input->getValue("measurePoints"));
    para->setpropellerValues(gridPath + input->getValue("propellerValues"));
    para->setclockCycleForMP(StringUtil::toFloat(input->getValue("measureClockCycle")));
    para->settimestepForMP(StringUtil::toInt(input->getValue("measureTimestep")));
    para->setcpTop(gridPath + input->getValue("cpTop"));
    para->setcpBottom(gridPath + input->getValue("cpBottom"));
    para->setcpBottom2(gridPath + input->getValue("cpBottom2"));
    para->setConcentration(gridPath + input->getValue("Concentration"));
    //////////////////////////////////////////////////////////////////////////
    //Normals - Geometry
    para->setgeomBoundaryNormalX(gridPath + input->getValue("geomBoundaryNormalX"));
    para->setgeomBoundaryNormalY(gridPath + input->getValue("geomBoundaryNormalY"));
    para->setgeomBoundaryNormalZ(gridPath + input->getValue("geomBoundaryNormalZ"));
    //Normals - Inlet
    para->setInflowBoundaryNormalX(gridPath + input->getValue("inletBoundaryNormalX"));
    para->setInflowBoundaryNormalY(gridPath + input->getValue("inletBoundaryNormalY"));
    para->setInflowBoundaryNormalZ(gridPath + input->getValue("inletBoundaryNormalZ"));
    //Normals - Outlet
    para->setOutflowBoundaryNormalX(gridPath + input->getValue("outletBoundaryNormalX"));
    para->setOutflowBoundaryNormalY(gridPath + input->getValue("outletBoundaryNormalY"));
    para->setOutflowBoundaryNormalZ(gridPath + input->getValue("outletBoundaryNormalZ"));
    //////////////////////////////////////////////////////////////////////////
    //Forcing
    para->setForcing(StringUtil::toFloat(input->getValue("ForcingX")), StringUtil::toFloat(input->getValue("ForcingY")), StringUtil::toFloat(input->getValue("ForcingZ")));
    //////////////////////////////////////////////////////////////////////////
    //Particles
    para->setCalcParticles(StringUtil::toBool(input->getValue("calcParticles")));
    para->setParticleBasicLevel(StringUtil::toInt(input->getValue("baseLevel")));
    para->setParticleInitLevel(StringUtil::toInt(input->getValue("initLevel")));
    para->setNumberOfParticles(StringUtil::toInt(input->getValue("numberOfParticles")));
    para->setneighborWSB(gridPath + input->getValue("neighborWSB"));
    para->setStartXHotWall(StringUtil::toDouble(input->getValue("startXHotWall")));
    para->setEndXHotWall(StringUtil::toDouble(input->getValue("endXHotWall")));
    //////////////////////////////////////////////////////////////////////////
    //for Multi GPU
    if (para->getNumprocs() > 1)
    {
        ////////////////////////////////////////////////////////////////////////////
        ////1D domain decomposition
        //std::vector<std::string> sendProcNeighbors;
        //std::vector<std::string> recvProcNeighbors;
        //for (int i = 0; i<para->getNumprocs();i++)
        //{
        // sendProcNeighbors.push_back(gridPath + StringUtil::toString(i) + "s.dat");
        // recvProcNeighbors.push_back(gridPath + StringUtil::toString(i) + "r.dat");
        //}
        //para->setPossNeighborFiles(sendProcNeighbors, "send");
        //para->setPossNeighborFiles(recvProcNeighbors, "recv");
        //////////////////////////////////////////////////////////////////////////
        //3D domain decomposition
        std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
        std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
        for (int i = 0; i < para->getNumprocs(); i++)
        {
            sendProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xs.dat");
            sendProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Ys.dat");
            sendProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zs.dat");
            recvProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xr.dat");
            recvProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Yr.dat");
            recvProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zr.dat");
        }
        para->setPossNeighborFilesX(sendProcNeighborsX, "send");
        para->setPossNeighborFilesY(sendProcNeighborsY, "send");
        para->setPossNeighborFilesZ(sendProcNeighborsZ, "send");
        para->setPossNeighborFilesX(recvProcNeighborsX, "recv");
        para->setPossNeighborFilesY(recvProcNeighborsY, "recv");
        para->setPossNeighborFilesZ(recvProcNeighborsZ, "recv");
    }
    //////////////////////////////////////////////////////////////////////////
    //para->setkFull(             input->getValue( "kFull" ));
    //para->setgeoFull(           input->getValue( "geoFull" ));
    //para->setnoSlipBcPos(       input->getValue( "noSlipBcPos" ));
    //para->setnoSlipBcQs(          input->getValue( "noSlipBcQs" ));
    //para->setnoSlipBcValues(      input->getValue( "noSlipBcValues" ));
    //para->setnoSlipBcValue(     input->getValue( "noSlipBcValue" ));
    //para->setslipBcPos(         input->getValue( "slipBcPos" ));
    //para->setslipBcQs(          input->getValue( "slipBcQs" ));
    //para->setslipBcValue(       input->getValue( "slipBcValue" ));
    //para->setpressBcPos(        input->getValue( "pressBcPos" ));
    //para->setpressBcQs(           input->getValue( "pressBcQs" ));
    //para->setpressBcValues(       input->getValue( "pressBcValues" ));
    //para->setpressBcValue(      input->getValue( "pressBcValue" ));
    //para->setvelBcQs(             input->getValue( "velBcQs" ));
    //para->setvelBcValues(         input->getValue( "velBcValues" ));
    //para->setpropellerCylinder( input->getValue( "propellerCylinder" ));
    //para->setpropellerQs(		 input->getValue( "propellerQs"      ));
    //para->setwallBcQs(            input->getValue( "wallBcQs"         ));
    //para->setwallBcValues(        input->getValue( "wallBcValues"     ));
    //para->setperiodicBcQs(        input->getValue( "periodicBcQs"     ));
    //para->setperiodicBcValues(    input->getValue( "periodicBcValues" ));
    //cout << "Try this: " << para->getgeomBoundaryBcValues() << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Restart
    para->setTimeDoCheckPoint(StringUtil::toInt(input->getValue("TimeDoCheckPoint")));
    para->setTimeDoRestart(StringUtil::toInt(input->getValue("TimeDoRestart")));
    para->setDoCheckPoint(StringUtil::toBool(input->getValue("DoCheckPoint")));
    para->setDoRestart(StringUtil::toBool(input->getValue("DoRestart")));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    para->setMaxLevel(StringUtil::toInt(input->getValue("NOGL")));
    para->setGridX(StringUtil::toIntVector(input->getValue("GridX")));                           
    para->setGridY(StringUtil::toIntVector(input->getValue("GridY")));                           
    para->setGridZ(StringUtil::toIntVector(input->getValue("GridZ")));                  
    para->setDistX(StringUtil::toIntVector(input->getValue("DistX")));                  
    para->setDistY(StringUtil::toIntVector(input->getValue("DistY")));                  
    para->setDistZ(StringUtil::toIntVector(input->getValue("DistZ")));                

    para->setNeedInterface(std::vector<bool>{true, true, true, true, true, true});
}



void multipleLevel(const std::string& configPath)
{
    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
    SPtr<Parameter> para = Parameter::make();
    SPtr<GridProvider> gridGenerator;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool useGridGenerator = true;

    if(useGridGenerator){
        
        //const uint generatePart = 1;
        const uint generatePart = Communicator::getInstanz()->getPID();
            
        real dx = 1.0 / 20.0;
        real vx = 0.05;

        auto triangularMesh = std::make_shared<TriangularMesh>("F:/Work/Computations/gridGenerator/stl/ShpereNotOptimal.stl");
        //auto triangularMesh = std::make_shared<TriangularMesh>("stl/ShpereNotOptimal.lnx.stl");

        // all
        //gridBuilder->addCoarseGrid(-2, -2, -2,  
        //                            4,  2,  2, dx);

        real overlap = 10.0 * dx;

        gridBuilder->addCoarseGrid(-2.0, -2.0, -2.0,  
                                    4.0,  2.0,  2.0, dx);


        gridBuilder->setNumberOfLayers(10,8);
        gridBuilder->addGrid(triangularMesh, 1);

        gridBuilder->addGeometry(triangularMesh);

        gridBuilder->setPeriodicBoundaryCondition(false, false, false);

        gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

        //////////////////////////////////////////////////////////////////////////

        gridBuilder->setVelocityBoundaryCondition(SideType::PY, vx , 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MY, vx , 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

        gridBuilder->setVelocityBoundaryCondition(SideType::MX, vx, 0.0, 0.0);
        gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);

        gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
        bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
        bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

        //////////////////////////////////////////////////////////////////////////
        gridBuilder->writeGridsToVtk("F:/Work/Computations/gridGenerator/grid/Test_");
        //gridBuilder->writeArrows    ("F:/Work/Computations/gridGenerator/grid/Test_Arrow");

        //SimulationFileWriter::write("F:/Work/Computations/gridGenerator/grid/", gridBuilder, FILEFORMAT::ASCII);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
        if(false)
        {

            auto getParentIndex = [&] (uint index, uint level) -> uint
            {
                SPtr<Grid> grid = gridBuilder->getGrid( level );

                if( level != 0 )
                {
                    real x, y, z;
                    grid->transIndexToCoords(index, x, y, z);

                    SPtr<Grid> coarseGrid = gridBuilder->getGrid(level - 1);

                    for (const auto dir : DistributionHelper::getDistribution27())
                    {
                        if (std::abs(dir[0]) < 0.5 || std::abs(dir[1]) < 0.5 || std::abs(dir[2]) < 0.5) continue;

                        real coarseX = x + dir[0] * 0.5 * grid->getDelta();
                        real coarseY = y + dir[1] * 0.5 * grid->getDelta();
                        real coarseZ = z + dir[2] * 0.5 * grid->getDelta();

                        // check if close enough to coarse grid coordinates
                        if( 0.01 * grid->getDelta() < std::abs(         (coarseGrid->getStartX() - coarseX) / grid->getDelta() 
                                                                - lround( (coarseGrid->getStartX() - coarseX) / grid->getDelta() ) ) ) continue;
                        if( 0.01 * grid->getDelta() < std::abs(         (coarseGrid->getStartY() - coarseY) / grid->getDelta() 
                                                                - lround( (coarseGrid->getStartY() - coarseY) / grid->getDelta() ) ) ) continue;
                        if( 0.01 * grid->getDelta() < std::abs(         (coarseGrid->getStartZ() - coarseZ) / grid->getDelta() 
                                                                - lround( (coarseGrid->getStartZ() - coarseZ) / grid->getDelta() ) ) ) continue;

                        uint parentIndex = coarseGrid->transCoordToIndex( coarseX, coarseY, coarseZ);

                        return parentIndex;
                    }
                }

                return INVALID_INDEX;
            };


            std::vector<idx_t> xadj;
            std::vector<idx_t> adjncy;

            std::vector<idx_t> vwgt;
            std::vector<idx_t> adjwgt;

            idx_t vertexCounter = 0;
            uint edgeCounter = 0;

            std::cout << "Checkpoint 1:" << std::endl;

            std::vector< std::vector<idx_t> > vertexIndex( gridBuilder->getNumberOfLevels() );

            std::vector< uint > startVerticesPerLevel;;

            for( uint level = 0; level < gridBuilder->getNumberOfLevels(); level++ )
            {
                SPtr<Grid> grid = gridBuilder->getGrid( level );

                vertexIndex[level].resize( grid->getSize() );

                startVerticesPerLevel.push_back(vertexCounter);

                for (uint index = 0; index < grid->getSize(); index++)
                {
                    if (grid->getSparseIndex(index) == INVALID_INDEX)
                    {
                        vertexIndex[level][index] = INVALID_INDEX;
                        continue;
                    }

                    uint parentIndex = getParentIndex(index, level);

                    if( parentIndex != INVALID_INDEX )
                    {
                        SPtr<Grid> coarseGrid = gridBuilder->getGrid(level - 1);

                        if( coarseGrid->getFieldEntry(parentIndex) == FLUID_CFC ||
                            coarseGrid->getFieldEntry(parentIndex) == FLUID_FCC ||
                            coarseGrid->getFieldEntry(parentIndex) == STOPPER_COARSE_UNDER_FINE )
                        {
                            //vertexIndex[level][index] = INVALID_INDEX;
                            vertexIndex[level][index] = vertexIndex[level - 1][parentIndex];
                            continue;
                        }
                    }

                    vertexIndex[level][index] = vertexCounter;

                    vwgt.push_back( std::pow(2, level) );
                    //vwgt.push_back( std::pow(2, 2*level) );
                    vertexCounter++;
                }

            }

            //////////////////////////////////////////////////////////////////////////
            //for( uint level = 0; level < gridBuilder->getNumberOfLevels(); level++ )
            //{
            //    SPtr<Grid> grid = gridBuilder->getGrid( level );

            //    for (uint index = 0; index < grid->getSize(); index++)
            //    {
            //        grid->setFieldEntry(index, vertexIndex[level][index] >= startVerticesPerLevel[level] && vertexIndex[level][index] != INVALID_INDEX);
            //    }
            //}

            //gridBuilder->writeGridsToVtk("F:/Work/Computations/gridGenerator/grid/VertexIndex_");

            //return;
            //////////////////////////////////////////////////////////////////////////


            std::cout << "Checkpoint 2:" << std::endl;
                
            for( uint level = 0; level < gridBuilder->getNumberOfLevels(); level++ )
            {
                SPtr<Grid> grid = gridBuilder->getGrid( level );

                for (uint index = 0; index < grid->getSize(); index++)
                {
                    //if (grid->getSparseIndex(index) == INVALID_INDEX) continue;

                    if( vertexIndex[level][index] == INVALID_INDEX ) continue;

                    if( vertexIndex[level][index] < startVerticesPerLevel[level] ) continue;

                    xadj.push_back(edgeCounter);

                    real x, y, z;
                    grid->transIndexToCoords(index, x, y, z);

                    for (const auto dir : DistributionHelper::getDistribution27())
                    {
                        const uint neighborIndex = grid->transCoordToIndex(x + dir[0] * grid->getDelta(), 
                                                                            y + dir[1] * grid->getDelta(), 
                                                                            z + dir[2] * grid->getDelta());

                        if (neighborIndex == INVALID_INDEX) continue;

                        if (neighborIndex == index) continue;

                        if( vertexIndex[level][neighborIndex] == INVALID_INDEX ) continue;

                        adjncy.push_back( vertexIndex[level][neighborIndex] );
                        adjwgt.push_back( std::pow(2, level) );

                        edgeCounter++;
                    }

                    if( grid->getFieldEntry(index) == FLUID_CFC ||
                        grid->getFieldEntry(index) == FLUID_FCC ||
                        grid->getFieldEntry(index) == STOPPER_COARSE_UNDER_FINE )

                    {
                        SPtr<Grid> fineGrid = gridBuilder->getGrid(level + 1);

                        for (const auto dir : DistributionHelper::getDistribution27())
                        {
                            if (std::abs(dir[0]) < 0.5 || std::abs(dir[1]) < 0.5 || std::abs(dir[2]) < 0.5) continue;

                            real fineX = x + dir[0] * 0.25 * grid->getDelta();
                            real fineY = y + dir[1] * 0.25 * grid->getDelta();
                            real fineZ = z + dir[2] * 0.25 * grid->getDelta();

                            uint childIndex = fineGrid->transCoordToIndex(fineX, fineY, fineZ);

                            if( fineGrid->getFieldEntry(childIndex) == INVALID_INDEX ) continue;
                            if( vertexIndex[level + 1][childIndex]  == INVALID_INDEX ) continue;

                            for (const auto dir : DistributionHelper::getDistribution27())
                            {
                                const uint neighborIndex = fineGrid->transCoordToIndex( fineX + dir[0] * fineGrid->getDelta(), 
                                                                                        fineY + dir[1] * fineGrid->getDelta(), 
                                                                                        fineZ + dir[2] * fineGrid->getDelta() );

                                if(neighborIndex == INVALID_INDEX) continue;

                                if (neighborIndex == childIndex) continue;

                                if( vertexIndex[level + 1][neighborIndex] == INVALID_INDEX ) continue;

                                adjncy.push_back( vertexIndex[level + 1][neighborIndex] );
                                adjwgt.push_back( std::pow(2, level) );

                                edgeCounter++;
                            }
                        }
                    }
                }
            }

            xadj.push_back( edgeCounter );

            std::cout << "Checkpoint 3:" << std::endl;
                
            idx_t nWeights  = 1;
            idx_t nParts    = 4;
            idx_t objval    = 0;

            std::vector<idx_t> part( vertexCounter );
                
            std::cout << vertexCounter << std::endl;
            std::cout << edgeCounter << std::endl;
            std::cout << xadj.size()  << std::endl;
            std::cout << adjncy.size() << std::endl;

            //int ret = METIS_PartGraphRecursive(&vertexCounter, &nWeights, xadj.data(), adjncy.data(),
            // 				                   vwgt.data(), NULL, adjwgt.data(), &nParts, 
            //                                   NULL, NULL, NULL, &objval, part.data());

            int ret = METIS_PartGraphKway(&vertexCounter, &nWeights, xadj.data(), adjncy.data(),
                 				            vwgt.data(), NULL, NULL/*adjwgt.data()*/, &nParts, 
                                            NULL, NULL, NULL, &objval, part.data());

            std::cout << "objval:" << objval << std::endl;

            std::cout << "Checkpoint 4:" << std::endl;

            //uint partCounter = 0;
                
            for( uint level = 0; level < gridBuilder->getNumberOfLevels(); level++ )
            {
                SPtr<Grid> grid = gridBuilder->getGrid( level );

                for (uint index = 0; index < grid->getSize(); index++)
                {
                    if (grid->getSparseIndex(index) == INVALID_INDEX) continue;

                    grid->setFieldEntry(index, part[vertexIndex[level][index]]);

                    //partCounter++;
                }
            }

            std::cout << "Checkpoint 5:" << std::endl;

            gridBuilder->writeGridsToVtk("F:/Work/Computations/gridGenerator/grid/Partition_");

        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
        {

            for( int level = gridBuilder->getNumberOfLevels()-1; level >= 0 ; level-- )
            {
                std::vector< std::vector<idx_t> > vertexIndex( gridBuilder->getNumberOfLevels() );

                std::vector<idx_t> xadj;
                std::vector<idx_t> adjncy;

                std::vector<idx_t> vwgt;
                std::vector<idx_t> adjwgt;

                idx_t vertexCounter = 0;
                uint edgeCounter = 0;

                SPtr<Grid> grid = gridBuilder->getGrid( level );

                vertexIndex[level].resize( grid->getSize() );

                for (uint index = 0; index < grid->getSize(); index++)
                {
                    if (grid->getSparseIndex(index) == INVALID_INDEX)
                    {
                        vertexIndex[level][index] = INVALID_INDEX;
                        continue;
                    }

                    vertexIndex[level][index] = vertexCounter;

                    vwgt.push_back( std::pow(2, level) );
                    //vwgt.push_back( std::pow(2, 2*level) );
                    vertexCounter++;
                }

                for (uint index = 0; index < grid->getSize(); index++)
                {
                    //if (grid->getSparseIndex(index) == INVALID_INDEX) continue;

                    if( vertexIndex[level][index] == INVALID_INDEX ) continue;

                    xadj.push_back(edgeCounter);

                    real x, y, z;
                    grid->transIndexToCoords(index, x, y, z);

                    for (const auto dir : DistributionHelper::getDistribution27())
                    {
                        const uint neighborIndex = grid->transCoordToIndex(x + dir[0] * grid->getDelta(), 
                                                                            y + dir[1] * grid->getDelta(), 
                                                                            z + dir[2] * grid->getDelta());

                        if (neighborIndex == INVALID_INDEX) continue;

                        if (neighborIndex == index) continue;

                        if( vertexIndex[level][neighborIndex] == INVALID_INDEX ) continue;

                        adjncy.push_back( vertexIndex[level][neighborIndex] );
                        adjwgt.push_back( std::pow(2, level) );

                        edgeCounter++;
                    }
                }

                xadj.push_back( edgeCounter );

                std::cout << "Checkpoint 3:" << std::endl;
                
                idx_t nWeights  = 1;
                idx_t nParts    = 4;
                idx_t objval    = 0;

                std::vector<idx_t> part( vertexCounter );
                
                std::cout << vertexCounter << std::endl;
                std::cout << edgeCounter << std::endl;
                std::cout << xadj.size()  << std::endl;
                std::cout << adjncy.size() << std::endl;

                int ret = METIS_PartGraphRecursive(&vertexCounter, &nWeights, xadj.data(), adjncy.data(),
                     				                NULL/*vwgt.data()*/, NULL, NULL/*adjwgt.data()*/, &nParts, 
                                                    NULL, NULL, NULL, &objval, part.data());

                //int ret = METIS_PartGraphKway(&vertexCounter, &nWeights, xadj.data(), adjncy.data(),
                 		//	                  NULL/*vwgt.data()*/, NULL, NULL/*adjwgt.data()*/, &nParts, 
                //                              NULL, NULL, NULL, &objval, part.data());

                std::cout << "objval:" << objval << std::endl;

                std::cout << "Checkpoint 4:" << std::endl;

                for (uint index = 0; index < grid->getSize(); index++)
                {
                    if (vertexIndex[level][index] == INVALID_INDEX) continue;

                    if( grid->getFieldEntry(index) == FLUID_CFC ||
                        grid->getFieldEntry(index) == FLUID_FCC ||
                        grid->getFieldEntry(index) == STOPPER_COARSE_UNDER_FINE )
                    {
                        SPtr<Grid> fineGrid = gridBuilder->getGrid(level+1);
                            
                        real x, y, z;
                        grid->transIndexToCoords(index, x, y, z);

                        for (const auto dir : DistributionHelper::getDistribution27())
                        {
                            if (std::abs(dir[0]) < 0.5 || std::abs(dir[1]) < 0.5 || std::abs(dir[2]) < 0.5) continue;

                            real fineX = x + dir[0] * 0.25 * grid->getDelta();
                            real fineY = y + dir[1] * 0.25 * grid->getDelta();
                            real fineZ = z + dir[2] * 0.25 * grid->getDelta();

                            uint childIndex = fineGrid->transCoordToIndex(fineX, fineY, fineZ);

                            if( childIndex == INVALID_INDEX ) continue;

                            fineGrid->setFieldEntry(childIndex, part[vertexIndex[level][index]]);
                            //fineGrid->setFieldEntry(childIndex, grid->getFieldEntry(index));
                        }
                    }

                    grid->setFieldEntry(index, part[vertexIndex[level][index]]);
                }
            }

            gridBuilder->writeGridsToVtk("F:/Work/Computations/gridGenerator/grid/Partition_");

        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        return;

        gridGenerator = GridGenerator::make(gridBuilder, para);
    }
    else
    {
        gridGenerator = GridReader::make(FileFormat::BINARY, para);
        //gridGenerator = GridReader::make(FileFormat::ASCII, para);
    }

    logFile.close();

    //return;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::ifstream stream;
    stream.open(configPath.c_str(), std::ios::in);
    if (stream.fail())
        throw std::runtime_error("can not open config file!");

    UPtr<input::Input> input = input::Input::makeInput(stream, "config");

    setParameters(para, input);

    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    sim.init(para, gridGenerator, fileWriter);
    sim.run();
}


int main( int argc, char* argv[])
{
     MPI_Init(&argc, &argv);
    std::string str, str2; 
    if ( argv != NULL )
    {
        str = static_cast<std::string>(argv[0]);
        if (argc > 1)
        {
            str2 = static_cast<std::string>(argv[1]);
            try
            {
                multipleLevel(str2);
            }
            catch (const std::exception& e)
            {
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (...)
            {
                std::cout << "unknown exeption" << std::endl;
            }
        }
        else
        {
            try
            {
                multipleLevel("F:/Work/Computations/gridGenerator/inp/configTest.txt");
            }
            catch (const std::exception& e)
            {
                std::cout << e.what() << std::flush;
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (const std::bad_alloc e)
            {
                std::cout << e.what() << std::flush;
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (...)
            {
                std::cout << "unknown exeption" << std::endl;
            }

            std::cout << "\nConfiguration file must be set!: lbmgm <config file>" << std::endl << std::flush;
            //MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }


   /*
   MPE_Init_log() & MPE_Finish_log() are NOT needed when
   liblmpe.a is linked with this program.  In that case,
   MPI_Init() would have called MPE_Init_log() already.
   */
#if defined( MPI_LOGGING )
   MPE_Init_log();
#endif

#if defined( MPI_LOGGING )
   if ( argv != NULL )
      MPE_Finish_log( argv[0] );
   if ( str != "" )
      MPE_Finish_log( str.c_str() );
   else
      MPE_Finish_log( "TestLog" );
#endif

   MPI_Finalize();
   return 0;
}
