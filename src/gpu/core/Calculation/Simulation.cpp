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
//! \author Martin Schoenherr
//=======================================================================================
#include "Simulation.h"

#include <vector>

#include <helper_timer.h>

#include "GridScaling/GridScalingFactory.h"
#include "Calculation/Calculation.h"
#include "Communication/ExchangeData27.h"
#include "Parameter/Parameter.h"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/EdgeNodeFinder.h"

#include "basics/utilities/UbFileOutputASCII.h"
//////////////////////////////////////////////////////////////////////////
#include "Output/MeasurePointWriter.hpp"
#include "Output/AnalysisData.hpp"
#include "Output/InterfaceDebugWriter.hpp"
#include "Output/EdgeNodeDebugWriter.hpp"
#include "Output/NeighborDebugWriter.hpp"
#include "Output/VeloASCIIWriter.hpp"
//////////////////////////////////////////////////////////////////////////
#include "StringUtilities/StringUtil.h"
//////////////////////////////////////////////////////////////////////////
#include "PreProcessor/InitLattice.h"
#include "PreProcessor/ReaderMeasurePoints.h"
//////////////////////////////////////////////////////////////////////////
#include "PostProcessor/Calc2ndMoments.h"
#include "PostProcessor/CalcMean.h"
#include "PostProcessor/CalcTurbulenceIntensity.h"
#include "PostProcessor/Concentration.cuh"
#include "PostProcessor/Cp.h"
#include "PostProcessor/DragLift.h"
#include "PostProcessor/EnstrophyAnalyzer.h"
#include "PostProcessor/ForceCalculations.h"
#include "PostProcessor/KineticEnergyAnalyzer.h"
#include "PostProcessor/MacroscopicQuantities.cuh"
#include "PostProcessor/PlaneCalculations.h"
#include "PostProcessor/TurbulenceIntensity.h"
#include "UpdateGrid27.h"
//////////////////////////////////////////////////////////////////////////
#include "Output/FileWriter.h"
#include "Output/DistributionDebugWriter.h"
//////////////////////////////////////////////////////////////////////////
#include "Restart/RestartObject.h"
//////////////////////////////////////////////////////////////////////////
#include "DataStructureInitializer/GridProvider.h"
#include "Output/DataWriter.h"
#include "Kernel/KernelFactory/KernelFactory.h"
#include "PreProcessor/PreProcessorFactory/PreProcessorFactory.h"
#include "PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"
#include "Kernel/KernelFactory/KernelFactoryImp.h"
#include "Kernel/Kernel.h"
#include "TurbulenceModels/TurbulenceModelFactory.h"

#include <cuda_helper/DeviceInfo.h>

#include <logger/Logger.h>

#include <parallel/Communicator.h>

std::string getFileName(const std::string& fname, int step, int myID)
{
    return std::string(fname + "_Restart_" + UbSystem::toString(myID) + "_" +  UbSystem::toString(step));
}

Simulation::Simulation(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> memoryManager,
                       vf::parallel::Communicator &communicator, GridProvider &gridProvider, BoundaryConditionFactory* bcFactory, GridScalingFactory* scalingFactory)
    : para(para), cudaMemoryManager(memoryManager), communicator(communicator), kernelFactory(std::make_unique<KernelFactoryImp>()),
      preProcessorFactory(std::make_shared<PreProcessorFactoryImp>()), dataWriter(std::make_unique<FileWriter>())
{
    this->tmFactory = SPtr<TurbulenceModelFactory>( new TurbulenceModelFactory(para) );
    init(gridProvider, bcFactory, tmFactory, scalingFactory);
}

Simulation::Simulation(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> memoryManager,
                       vf::parallel::Communicator &communicator, GridProvider &gridProvider, BoundaryConditionFactory* bcFactory, SPtr<TurbulenceModelFactory> tmFactory, GridScalingFactory* scalingFactory)
    : para(para), cudaMemoryManager(memoryManager), communicator(communicator), kernelFactory(std::make_unique<KernelFactoryImp>()),
      preProcessorFactory(std::make_shared<PreProcessorFactoryImp>()), dataWriter(std::make_unique<FileWriter>())
{
    init(gridProvider, bcFactory, tmFactory, scalingFactory);
}

void Simulation::init(GridProvider &gridProvider, BoundaryConditionFactory *bcFactory, SPtr<TurbulenceModelFactory> tmFactory, GridScalingFactory *scalingFactory)
{
    gridProvider.initalGridInformations();

    vf::cuda::verifyAndSetDevice(communicator.mapCudaDevicesOnHosts(para->getDevices(), para->getMaxDev()));

    para->initLBMSimulationParameter();

    gridProvider.allocAndCopyForcing();
    gridProvider.allocAndCopyQuadricLimiters();
        
    gridProvider.setDimensions();
    gridProvider.setBoundingBox();

    para->setRe(para->getVelocity() * (real)1.0 / para->getViscosity());
    if (para->getDoRestart())
        para->setStartTurn(para->getTimeDoRestart());
    else
        para->setStartTurn((unsigned int)0);

    restart_object = std::make_shared<ASCIIRestartObject>();

    //////////////////////////////////////////////////////////////////////////
    VF_LOG_INFO("LB_Modell:       D3Q{}", para->getD3Qxx());
    VF_LOG_INFO("Re:              {}", para->getRe());
    VF_LOG_INFO("vis_ratio:       {}", para->getViscosityRatio());
    VF_LOG_INFO("u0_ratio:        {}", para->getVelocityRatio());
    VF_LOG_INFO("delta_rho:       {}", para->getDensityRatio());
    VF_LOG_INFO("QuadricLimiters: {}, \t{}, \t{}", para->getQuadricLimitersHost()[0],
                para->getQuadricLimitersHost()[1], para->getQuadricLimitersHost()[2]);
    //////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    cudaMemoryManager->setMemsizeGPU(0, true);
    //////////////////////////////////////////////////////////////////////////
    allocNeighborsOffsetsScalesAndBoundaries(gridProvider);

    //! Get tagged fluid nodes with corresponding value for CollisionTemplate from interactors
    for (SPtr<PreCollisionInteractor> actuator : para->getActuators()) {
        actuator->init(para.get(), &gridProvider, cudaMemoryManager.get());
        actuator->getTaggedFluidNodes( para.get(), &gridProvider );
    }

    for (SPtr<PreCollisionInteractor> probe : para->getProbes()) {
        probe->init(para.get(), &gridProvider, cudaMemoryManager.get());
        probe->getTaggedFluidNodes( para.get(), &gridProvider );
    }

    //////////////////////////////////////////////////////////////////////////
    // CUDA streams
    if (para->getUseStreams()) {
        para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::SubDomainBorder);
        para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::Bulk);
        para->getStreamManager()->createCudaEvents();
    }
    //////////////////////////////////////////////////////////////////////////
    
    if (para->getKernelNeedsFluidNodeIndicesToRun())
    {
        gridProvider.sortFluidNodeTags();
        gridProvider.allocArrays_taggedFluidNodes();
    }
    //////////////////////////////////////////////////////////////////////////
    // Kernel init
    //////////////////////////////////////////////////////////////////////////
    VF_LOG_INFO("make Kernels");
    kernels = kernelFactory->makeKernels(para);

    if (para->getDiffOn()) {
        VF_LOG_INFO("make AD Kernels");
        adKernels = kernelFactory->makeAdvDifKernels(para);
        std::vector<PreProcessorType> preProADTypes = adKernels.at(0)->getPreProcessorTypes();
        preProcessorAD = preProcessorFactory->makePreProcessor(preProADTypes, para);
    }

    //////////////////////////////////////////////////////////////////////////
    // PreProcessor init
    //////////////////////////////////////////////////////////////////////////
    VF_LOG_INFO("make Preprocessors");
    std::vector<PreProcessorType> preProTypes = kernels.at(0)->getPreProcessorTypes();
    preProcessor = preProcessorFactory->makePreProcessor(preProTypes, para);

    //////////////////////////////////////////////////////////////////////////
    // Allocate Memory for Drag Lift Calculation
    //////////////////////////////////////////////////////////////////////////
    if (para->getCalcDragLift())
        allocDragLift(para.get(), cudaMemoryManager.get());

    //////////////////////////////////////////////////////////////////////////
    // Allocate Memory for Plane Conc Calculation
    //////////////////////////////////////////////////////////////////////////
    // if (para->getDiffOn()) allocPlaneConc(para.get(), cudaMemoryManager.get());

    //////////////////////////////////////////////////////////////////////////
    // Mean
    //////////////////////////////////////////////////////////////////////////
    if (para->getCalcMean()) {
        VF_LOG_INFO("alloc Calculation for Mean Values");
        if (para->getDiffOn())
            allocMeanAD(para.get(), cudaMemoryManager.get());
        else
            allocMean(para.get(), cudaMemoryManager.get());
    }

    //////////////////////////////////////////////////////////////////////////
    // Turbulence Intensity
    //////////////////////////////////////////////////////////////////////////
    if (para->getCalcTurbulenceIntensity()) {
        VF_LOG_INFO("alloc arrays for calculating Turbulence Intensity");
        allocTurbulenceIntensity(para.get(), cudaMemoryManager.get());
    }

    //////////////////////////////////////////////////////////////////////////
    // allocate memory and initialize 2nd, 3rd and higher order moments
    //////////////////////////////////////////////////////////////////////////
    if (para->getCalc2ndOrderMoments()) {
        alloc2ndMoments(para.get(), cudaMemoryManager.get());
        init2ndMoments(para.get());
    }
    if (para->getCalc3rdOrderMoments()) {
        alloc3rdMoments(para.get(), cudaMemoryManager.get());
        init3rdMoments(para.get());
    }
    if (para->getCalcHighOrderMoments()) {
        allocHigherOrderMoments(para.get(), cudaMemoryManager.get());
        initHigherOrderMoments(para.get());
    }

    //////////////////////////////////////////////////////////////////////////
    // MeasurePoints
    //////////////////////////////////////////////////////////////////////////
    if (para->getUseMeasurePoints()) {
        VF_LOG_INFO("read measure points");
        ReaderMeasurePoints::readMeasurePoints(para.get(), cudaMemoryManager.get());
    }

    //////////////////////////////////////////////////////////////////////////
    // enSightGold
    //////////////////////////////////////////////////////////////////////////
    // excludeGridInterfaceNodesForMirror(para, 7);
    ////VF_LOG_INFO("print case file...");
    // printCaseFile(para);
    ////VF_LOG_INFO("print geo file...");
    // printGeoFile(para, true);  //true for binary
    ////VF_LOG_INFO("done.");

    //////////////////////////////////////////////////////////////////////////
    // Forcing
    //////////////////////////////////////////////////////////////////////////
    ////allocVeloForForcing(para);
    // VF_LOG_INFO("new object forceCalulator");
    // forceCalculator = std::make_shared<ForceCalculations>(para.get());

    //////////////////////////////////////////////////////////////////////////
    // VF_LOG_INFO("define the Grid...");
    // defineGrid(para, communicator);
    ////allocateMemory();
    // VF_LOG_INFO("done.");

    VF_LOG_INFO("init lattice...");
    initLattice(para, preProcessor, preProcessorAD, cudaMemoryManager);
    VF_LOG_INFO("done");

    // VF_LOG_INFO("set geo for Q...\n");
    // setGeoForQ();

    // if (maxlevel>1)
    //{
    // VF_LOG_INFO("find Qs...");
    // findQ27(para);
    // VF_LOG_INFO("done.");
    //}

    // if (para->getDiffOn()==true)
    //{
    //   VF_LOG_INFO("define TempBC...");
    //   findTempSim(para);

    //   VF_LOG_INFO("define TempVelBC...");
    //   findTempVelSim(para);

    //   VF_LOG_INFO("define TempPressBC...");
    //   findTempPressSim(para);
    //   VF_LOG_INFO("done.");
    //}

    // VF_LOG_INFO("find Qs-BC...");
    // findBC27(para);

    // VF_LOG_INFO("find Press-BC...");
    // findPressQShip(para);
    // VF_LOG_INFO("done.");

    //////////////////////////////////////////////////////////////////////////
    // find indices of corner nodes for multiGPU communication
    //////////////////////////////////////////////////////////////////////////
    if (para->getDevices().size() > 2) {
        VF_LOG_INFO("Find indices of edge nodes for multiGPU communication");
        vf::gpu::findEdgeNodesCommMultiGPU(*para);
    }
    //////////////////////////////////////////////////////////////////////////
    // Memory alloc for CheckPoint / Restart
    //////////////////////////////////////////////////////////////////////////
    if (para->getDoCheckPoint() || para->getDoRestart()) {
        VF_LOG_INFO("Alloc Memory for CheckPoint / Restart");
        for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
            cudaMemoryManager->cudaAllocFsForCheckPointAndRestart(lev);
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // Restart
    //////////////////////////////////////////////////////////////////////////
    if (para->getDoRestart()) {
        VF_LOG_INFO("Restart...\n...get the Object...");

        const auto name = getFileName(para->getFName(), para->getTimeDoRestart(), para->getMyProcessID());
        restart_object->deserialize(name, para);

        VF_LOG_INFO("...copy Memory for Restart...");
        for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
            //////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaCopyFsForRestart(lev);
            //////////////////////////////////////////////////////////////////////////
            // macroscopic values
            calculateMacroscopicQuantities(para->getParD(lev)->velocityX, para->getParD(lev)->velocityY, para->getParD(lev)->velocityZ,
                        para->getParD(lev)->rho, para->getParD(lev)->pressure, para->getParD(lev)->typeOfGridNode,
                        para->getParD(lev)->neighborX, para->getParD(lev)->neighborY,
                        para->getParD(lev)->neighborZ, para->getParD(lev)->numberOfNodes,
                        para->getParD(lev)->numberofthreads, para->getParD(lev)->distributions.f[0],
                        para->getParD(lev)->isEvenTimestep);
            getLastCudaError("Kernel calculateMacroscopicQuantities execution failed");
            //////////////////////////////////////////////////////////////////////////
            // test...should not work...and does not
            // para->getEvenOrOdd(lev)==false;
        }
        VF_LOG_INFO("done.");
    }

    //////////////////////////////////////////////////////////////////////////
    // Init UpdateGrid
    //////////////////////////////////////////////////////////////////////////
    this->updateGrid27 = std::make_unique<UpdateGrid27>(para, communicator, cudaMemoryManager, kernels, adKernels, bcFactory, tmFactory, scalingFactory);

    //////////////////////////////////////////////////////////////////////////
    // Write Initialized Files
    //////////////////////////////////////////////////////////////////////////
    VF_LOG_INFO("Write initialized Files ...");
    dataWriter->writeInit(para, cudaMemoryManager);
    VF_LOG_INFO("... done.");

    // VF_LOG_INFO("Write vtk files for debugging...");
    // NeighborDebugWriter::writeNeighborLinkLinesDebug(para.get());

    // InterfaceDebugWriter::writeInterfaceLinesDebugCF(para.get());
    // InterfaceDebugWriter::writeInterfaceLinesDebugFC(para.get());

    // writers for version with communication hiding
    //    if(para->getNumprocs() > 1 && para->getUseStreams()){
    //        InterfaceDebugWriter::writeInterfaceFCC_Send(para.get());
    //        InterfaceDebugWriter::writeInterfaceCFC_Recv(para.get());
    //        InterfaceDebugWriter::writeSendNodesStream(para.get());
    //        InterfaceDebugWriter::writeRecvNodesStream(para.get());
    //        EdgeNodeDebugWriter::writeEdgeNodesXZ_Send(para);
    //        EdgeNodeDebugWriter::writeEdgeNodesXZ_Recv(para);
    //    }

#if DEBUG_FS
    VF_LOG_INFO("Allocating host memory for DistributionWriter.");
    DistributionDebugWriter::allocateDistributionsOnHost(*cudaMemoryManager);
#endif

    // VF_LOG_INFO("...done");

    //////////////////////////////////////////////////////////////////////////
    VF_LOG_INFO("used Device Memory: {} MB", cudaMemoryManager->getMemsizeGPU() / 1000000.0);
    // std::cout << "Process " << communicator.getPID() <<": used device memory" << cudaMemoryManager->getMemsizeGPU() /
    // 1000000.0 << " MB\n" << std::endl;
    //////////////////////////////////////////////////////////////////////////
}

void Simulation::addKineticEnergyAnalyzer(uint tAnalyse)
{
    this->kineticEnergyAnalyzer = std::make_unique<KineticEnergyAnalyzer>(this->para, tAnalyse);
}

void Simulation::addEnstrophyAnalyzer(uint tAnalyse)
{
    this->enstrophyAnalyzer = std::make_unique<EnstrophyAnalyzer>(this->para, tAnalyse);
}

void Simulation::setDataWriter(std::shared_ptr<DataWriter> dataWriter_)
{
    this->dataWriter = dataWriter_;
}

void Simulation::setFactories(std::unique_ptr<KernelFactory> &&kernelFactory_,
               std::unique_ptr<PreProcessorFactory> &&preProcessorFactory_)
{
    this->kernelFactory = std::move(kernelFactory_);
    this->preProcessorFactory = std::move(preProcessorFactory_);
}

void Simulation::initTimers()
{
    previousTimestepForAveraging = para->getTimeCalcMedStart();
    previousTimestepForTurbulenceIntensityCalculation = 0;
    timestepForMeasuringPoints = 0;
    
    para->setStepEnsight(0);

    averageTimer.start();
}

void Simulation::allocNeighborsOffsetsScalesAndBoundaries(GridProvider &gridProvider)
{
    gridProvider.allocArrays_CoordNeighborGeo();
    gridProvider.allocArrays_OffsetScale();
    gridProvider.allocArrays_BoundaryValues(); // allocArrays_BoundaryValues() has to be called after allocArrays_OffsetScale() because of initCommunicationArraysForCommAfterFinetoCoarse()
    gridProvider.allocArrays_BoundaryQs();
}


void Simulation::run()
{
    this->initTimers();

    //////////////////////////////////////////////////////////////////////////
    // turning Ship
    real Pi = (real)3.14159265358979323846;
    real delta_x_F = (real)0.1;
    real delta_t_F = (real)((double)para->getVelocity() * (double)delta_x_F / (double)3.75);
    real delta_t_C = (real)(delta_t_F * pow(2., para->getMaxLevel()));
    real timesteps_C = (real)(12.5 / delta_t_C);
    real AngularVelocity = (real)(12.5 / timesteps_C * Pi / 180.);
    para->setAngularVelocity(AngularVelocity);
    for (int i = 0; i <= para->getMaxLevel(); i++) {
        para->getParD(i)->deltaPhi = (real)(para->getAngularVelocity() / (pow(2., i)));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Time loop
    ////////////////////////////////////////////////////////////////////////////////
    for(uint timestep=para->getTimestepStart();timestep<=para->getTimestepEnd();timestep++)
    {
        this->calculateTimestep(timestep);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // printDragLift(para);
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    if (para->getDiffOn() == true) printPlaneConc(para.get(), cudaMemoryManager.get());
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    ////for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    ////{
    ////    if (para->getParH(lev)->cpTop.size() > 0)
    ////    {
    ////        printCpTop(para, lev);
    ////    }
    ////}
    // for (int lev = 7; lev <= 8; lev++)
    //{
    //    printCpTop(para, lev);
    //}
    ////printCpTop(para);
    ////printCpBottom(para);
    ////printCpBottom2(para);
    ////////////////////////////////////////////////////////////////////////////////

    //  //////////////////////////////////////////////////////////////////////////
    //  //Copy Measure Values
    // for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    //{
    //    VF_LOG_INFO("Copy MeasurePoints at level = {}", lev);
    //    para->cudaCopyMeasurePointsToHost(lev);
    //    para->copyMeasurePointsArrayToVector(lev);
    //    VF_LOG_INFO("Write MeasurePoints at level = {}", lev);
    //    for(int j = 0; j < (int)para->getParH(lev)->MP.size(); j++)
    //    {
    //        MeasurePointWriter::writeMeasurePoints(para, lev, j, 0);
    //    }
    //}
    //  //////////////////////////////////////////////////////////////////////////
}

void Simulation::calculateTimestep(uint timestep)
{
    this->updateGrid27->updateGrid(0, timestep);
    ////////////////////////////////////////////////////////////////////////////////
    // run Analyzers for kinetic energy and enstrophy for TGV in 3D
    // these analyzers only work on level 0
    ////////////////////////////////////////////////////////////////////////////////
    if (this->kineticEnergyAnalyzer || this->enstrophyAnalyzer) {
        updateGrid27->exchangeData(0);
    }
    if( this->kineticEnergyAnalyzer ) this->kineticEnergyAnalyzer->run(timestep);
    if( this->enstrophyAnalyzer     ) this->enstrophyAnalyzer->run(timestep);
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //Calc Mean
    ////////////////////////////////////////////////////////////////////////////////
    if (para->getCalcMean() && ((int)timestep >= para->getTimeCalcMedStart()) && ((int)timestep <= para->getTimeCalcMedEnd()))
    {
        for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
        {
            //calculateMean(para->getParD(lev)->vx_SP_Med,
                  //      para->getParD(lev)->vy_SP_Med,
                  //      para->getParD(lev)->vz_SP_Med,
                  //      para->getParD(lev)->rho_SP_Med,
                  //      para->getParD(lev)->press_SP_Med,
                  //      para->getParD(lev)->geoSP,
                  //      para->getParD(lev)->neighborX_SP,
                  //      para->getParD(lev)->neighborY_SP,
                  //      para->getParD(lev)->neighborZ_SP,
                  //      para->getParD(lev)->size_Mat_SP,
                  //      para->getParD(lev)->numberofthreads,
                  //      para->getParD(lev)->d0SP.f[0],
                  //      para->getParD(lev)->evenOrOdd);
            //getLastCudaError("calculateMacroscopicQuantities execution failed");
            calculateMeanCompressible(para->getParD(lev)->vx_SP_Med,
                            para->getParD(lev)->vy_SP_Med,
                            para->getParD(lev)->vz_SP_Med,
                            para->getParD(lev)->rho_SP_Med,
                            para->getParD(lev)->press_SP_Med,
                            para->getParD(lev)->typeOfGridNode,
                            para->getParD(lev)->neighborX,
                            para->getParD(lev)->neighborY,
                            para->getParD(lev)->neighborZ,
                            para->getParD(lev)->numberOfNodes,
                            para->getParD(lev)->numberofthreads,
                            para->getParD(lev)->distributions.f[0],
                            para->getParD(lev)->isEvenTimestep);
            getLastCudaError("CalcMacMedCompSP27 execution failed");
        }
    }
    if (para->getCalcTurbulenceIntensity()) {
        for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
            CalcTurbulenceIntensityDevice(
                para->getParD(lev)->vxx,
                para->getParD(lev)->vyy,
                para->getParD(lev)->vzz,
                para->getParD(lev)->vxy,
                para->getParD(lev)->vxz,
                para->getParD(lev)->vyz,
                para->getParD(lev)->vx_mean,
                para->getParD(lev)->vy_mean,
                para->getParD(lev)->vz_mean,
                para->getParD(lev)->distributions.f[0],
                para->getParD(lev)->typeOfGridNode,
                para->getParD(lev)->neighborX,
                para->getParD(lev)->neighborY,
                para->getParD(lev)->neighborZ,
                para->getParD(lev)->numberOfNodes,
                para->getParD(lev)->isEvenTimestep,
                para->getParD(lev)->numberofthreads
            );
        }
    }
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // CheckPoint
    ////////////////////////////////////////////////////////////////////////////////
    if(para->getDoCheckPoint() && para->getTimeDoCheckPoint()>0 && timestep%para->getTimeDoCheckPoint()==0 && timestep>0 && !para->overWritingRestart(timestep))
    {
        averageTimer.end();
        //////////////////////////////////////////////////////////////////////////
        if( para->getDoCheckPoint() )
        {
            VF_LOG_INFO("Copy data for CheckPoint t = {}....", timestep);
            for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
            {
                cudaMemoryManager->cudaCopyFsForCheckPoint(lev);
            }
            VF_LOG_INFO("Write data for CheckPoint t = {}...", timestep);
            const auto name = getFileName(para->getFName(), timestep, para->getMyProcessID());
            restart_object->serialize(name, para);
            VF_LOG_INFO("done");
        }
        //////////////////////////////////////////////////////////////////////////
        averageTimer.start();
    }
    //////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //Measure Points
    ////////////////////////////////////////////////////////////////////////////////
    //set MP-Time
    if (para->getUseMeasurePoints())
    {
        if ((timestep%para->getTimestepForMP()) == 0)
        {
            unsigned int valuesPerClockCycle = (unsigned int)(para->getclockCycleForMP() / para->getTimestepForMP());
            for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
            {
                // VF_LOG_INFO("start level = {}", lev);
                calculateMeasurePoints(  para->getParD(lev)->VxMP,            para->getParD(lev)->VyMP,                 para->getParD(lev)->VzMP,
                                        para->getParD(lev)->RhoMP,           para->getParD(lev)->kMP,                  para->getParD(lev)->numberOfPointskMP,
                                        valuesPerClockCycle,                 timestepForMeasuringPoints,                                     para->getParD(lev)->typeOfGridNode,
                                        para->getParD(lev)->neighborX,       para->getParD(lev)->neighborY,            para->getParD(lev)->neighborZ,
                                        para->getParD(lev)->numberOfNodes,   para->getParD(lev)->distributions.f[0],   para->getParD(lev)->numberofthreads,
                                        para->getParD(lev)->isEvenTimestep);
            }
            timestepForMeasuringPoints++;
        }
        //Copy Measure Values
        if ((timestep % (unsigned int)para->getclockCycleForMP()) == 0)
        {
            for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
            {
                cudaMemoryManager->cudaCopyMeasurePointsToHost(lev);
                para->copyMeasurePointsArrayToVector(lev);
                VF_LOG_INFO("Write MeasurePoints at level = {} and timestep = {}", lev, timestep);
                for (int j = 0; j < (int)para->getParH(lev)->MP.size(); j++)
                {
                    MeasurePointWriter::writeMeasurePoints(para.get(), lev, j, timestep);
                }
                //MeasurePointWriter::calcAndWriteMeanAndFluctuations(para.get(), lev, t, para->getTStartOut());
            }
            timestepForMeasuringPoints = 0;
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    ////get concentration at the plane
    //////////////////////////////////////////////////////////////////////////////////
    if (para->getDiffOn() && para->getCalcPlaneConc())
    {
        PlaneConcThS27( para->getParD(0)->ConcPlaneIn,
                       para->getParD(0)->cpTopIndex,
                       para->getParD(0)->numberOfPointsCpTop,
                       para->getParD(0)->typeOfGridNode,
                       para->getParD(0)->neighborX,
                       para->getParD(0)->neighborY,
                       para->getParD(0)->neighborZ,
                       para->getParD(0)->numberOfNodes,
                       para->getParD(0)->numberofthreads,
                       para->getParD(0)->distributionsAD.f[0],
                       para->getParD(0)->isEvenTimestep);
        getLastCudaError("PlaneConcThS27 execution failed");
        PlaneConcThS27( para->getParD(0)->ConcPlaneOut1,
                        para->getParD(0)->cpBottomIndex,
                        para->getParD(0)->numberOfPointsCpBottom,
                        para->getParD(0)->typeOfGridNode,
                        para->getParD(0)->neighborX,
                        para->getParD(0)->neighborY,
                        para->getParD(0)->neighborZ,
                        para->getParD(0)->numberOfNodes,
                        para->getParD(0)->numberofthreads,
                        para->getParD(0)->distributionsAD.f[0],
                        para->getParD(0)->isEvenTimestep);
        getLastCudaError("PlaneConcThS27 execution failed");
        PlaneConcThS27( para->getParD(0)->ConcPlaneOut2,
                        para->getParD(0)->pressureBC.kN,
                        para->getParD(0)->pressureBC.numberOfBCnodes,
                        para->getParD(0)->typeOfGridNode,
                        para->getParD(0)->neighborX,
                        para->getParD(0)->neighborY,
                        para->getParD(0)->neighborZ,
                        para->getParD(0)->numberOfNodes,
                        para->getParD(0)->numberofthreads,
                        para->getParD(0)->distributionsAD.f[0],
                        para->getParD(0)->isEvenTimestep);
        getLastCudaError("PlaneConcThS27 execution failed");
        //////////////////////////////////////////////////////////////////////////////////
        ////Calculation of concentration at the plane
        //////////////////////////////////////////////////////////////////////////////////
        calcPlaneConc(para.get(), cudaMemoryManager.get(), 0);
    }
    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // File IO
    ////////////////////////////////////////////////////////////////////////////////
    if (para->getTimestepOut() > 0 && timestep % para->getTimestepOut() == 0 && timestep >= para->getTimestepStartOut()) {
        averageTimer.end();
        performanceOutput.print(averageTimer, timestep, para.get(), communicator);

        if (para->getPrintFiles()) {
            readAndWriteFiles(timestep);
        }
        averageTimer.start();
    }
}

void Simulation::readAndWriteFiles(uint timestep)
{
    VF_LOG_INFO("Write files t = {} ...", timestep);

    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //exchange data for valid post process
        updateGrid27->exchangeData(lev);

        // ////////////////////////////////////////////////////////////////////////
        // if (para->getD3Qxx()==19)
        // {
        //     CalcMac(para->getParD(lev)->vx,     para->getParD(lev)->vy,       para->getParD(lev)->vz,      para->getParD(lev)->rho,
        //             para->getParD(lev)->geo,    para->getParD(lev)->size_Mat, para->getParD(lev)->gridNX,  para->getParD(lev)->gridNY,
        //             para->getParD(lev)->gridNZ, para->getParD(lev)->d0.f[0],  para->getParD(lev)->evenOrOdd);
        // }
        // else if (para->getD3Qxx()==27)
        // {
        //    if (para->getCalcMean() && ((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()))
        //    {
        //         unsigned int tdiff = t - t_prev;
        //         calculateMacrosopicMean(para->getParD(lev)->vx_SP_Med,
        //                           para->getParD(lev)->vy_SP_Med,
        //                           para->getParD(lev)->vz_SP_Med,
        //                           para->getParD(lev)->rho_SP_Med,
        //                           para->getParD(lev)->press_SP_Med,
        //                           para->getParD(lev)->geoSP,
        //                           para->getParD(lev)->neighborX_SP,
        //                           para->getParD(lev)->neighborY_SP,
        //                           para->getParD(lev)->neighborZ_SP,
        //                           tdiff,
        //                           para->getParD(lev)->size_Mat_SP,
        //                           para->getParD(lev)->numberofthreads,
        //                           para->getParD(lev)->evenOrOdd);
        //         getLastCudaError("calculateMacrosopicMean execution failed");
        //    }
        //    calculateMacroscopicQuantities(para->getParD(lev)->vx_SP,
        //                        para->getParD(lev)->vy_SP,
        //                        para->getParD(lev)->vz_SP,
        //                        para->getParD(lev)->rho,
        //                        para->getParD(lev)->pressure,
        //                        para->getParD(lev)->geoSP,
        //                        para->getParD(lev)->neighborX_SP,
        //                        para->getParD(lev)->neighborY_SP,
        //                        para->getParD(lev)->neighborZ_SP,
        //                        para->getParD(lev)->size_Mat_SP,
        //                        para->getParD(lev)->numberofthreads,
        //                        para->getParD(lev)->d0SP.f[0],
        //                        para->getParD(lev)->evenOrOdd);
        //     getLastCudaError("calculateMacroscopicQuantities execution failed");
            calculateMacroscopicQuantitiesCompressible(para->getParD(lev)->velocityX,
                            para->getParD(lev)->velocityY,
                            para->getParD(lev)->velocityZ,
                            para->getParD(lev)->rho,
                            para->getParD(lev)->pressure,
                            para->getParD(lev)->typeOfGridNode,
                            para->getParD(lev)->neighborX,
                            para->getParD(lev)->neighborY,
                            para->getParD(lev)->neighborZ,
                            para->getParD(lev)->numberOfNodes,
                            para->getParD(lev)->numberofthreads,
                            para->getParD(lev)->distributions.f[0],
                            para->getParD(lev)->isEvenTimestep);
            getLastCudaError("calculateMacroscopicQuantities execution failed");

        cudaMemoryManager->cudaCopyPrint(lev);
        if (para->getCalcMean())
        {
            cudaMemoryManager->cudaCopyMeanPrint(lev);
        }
        //////////////////////////////////////////////////////////////////////////
        //TODO: implement flag to write ASCII data
        if (para->getWriteVeloASCIIfiles())
            VeloASCIIWriter::writeVelocitiesAsTXT(para.get(), lev, timestep);
        //////////////////////////////////////////////////////////////////////////
        if( this->kineticEnergyAnalyzer || this->enstrophyAnalyzer )
        {
            std::string fname = para->getFName() + "_ID_" + StringUtil::toString<int>(para->getMyProcessID()) + "_t_" + StringUtil::toString<int>(timestep);
            if (this->kineticEnergyAnalyzer) this->kineticEnergyAnalyzer->writeToFile(fname);
            if (this->enstrophyAnalyzer)     this->enstrophyAnalyzer->writeToFile(fname);
        }
        //////////////////////////////////////////////////////////////////////////
        if (para->getDiffOn())
        {
               CalcConcentration27(
                              para->getParD(lev)->numberofthreads,
                              para->getParD(lev)->concentration,
                              para->getParD(lev)->typeOfGridNode,
                              para->getParD(lev)->neighborX,
                              para->getParD(lev)->neighborY,
                              para->getParD(lev)->neighborZ,
                              para->getParD(lev)->numberOfNodes,
                              para->getParD(lev)->distributionsAD.f[0],
                              para->getParD(lev)->isEvenTimestep);
            cudaMemoryManager->cudaCopyConcentrationDeviceToHost(lev);
            //cudaMemoryCopy(para->getParH(lev)->Conc, para->getParD(lev)->Conc,  para->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost);
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////print cp
        //if ((para->getParH(lev)->cpTop.size() > 0) && (t > para->getTStartOut()))
        //{
           // printCpTopIntermediateStep(para, t, lev);
        //}
        ////////////////////////////////////////////////////////////////////////////////
        //MeasurePointWriter::writeSpacialAverageForXZSlices(para, lev, t);
        ////////////////////////////////////////////////////////////////////////////////
        //MeasurePointWriter::writeTestAcousticXY(para, lev, t);
        //MeasurePointWriter::writeTestAcousticYZ(para, lev, t);
        //MeasurePointWriter::writeTestAcousticXZ(para, lev, t);
        ////////////////////////////////////////////////////////////////////////
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////test print press mirror
    //if (t > para->getTStartOut())
    //{
    //    ////////////////////////////////////////////////////////////////////////////////
    //    //Level 7
    //    CalcCPtop27(para->getParD(7)->d0SP.f[0],
    //        para->getParD(7)->cpTopIndex,
    //        para->getParD(7)->numberOfPointsCpTop,
    //        para->getParD(7)->cpPressTop,
    //        para->getParD(7)->neighborX_SP,
    //        para->getParD(7)->neighborY_SP,
    //        para->getParD(7)->neighborZ_SP,
    //        para->getParD(7)->size_Mat_SP,
    //        para->getParD(7)->evenOrOdd,
    //        para->getParD(7)->numberofthreads);
    //    //////////////////////////////////////////////////////////////////////////////////
    //    calcPressForMirror(para, 7);
    //    ////////////////////////////////////////////////////////////////////////////////
    //    //Level 8
    //    CalcCPtop27(para->getParD(8)->d0SP.f[0],
    //        para->getParD(8)->cpTopIndex,
    //        para->getParD(8)->numberOfPointsCpTop,
    //        para->getParD(8)->cpPressTop,
    //        para->getParD(8)->neighborX_SP,
    //        para->getParD(8)->neighborY_SP,
    //        para->getParD(8)->neighborZ_SP,
    //        para->getParD(8)->size_Mat_SP,
    //        para->getParD(8)->evenOrOdd,
    //        para->getParD(8)->numberofthreads);
    //    //////////////////////////////////////////////////////////////////////////////////
    //    calcPressForMirror(para, 8);
    //    ////////////////////////////////////////////////////////////////////////////////
    //    //print press mirror
    //    printScalars(para, false);
    //    ////////////////////////////////////////////////////////////////////////////////
    //}
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //t_prev = t;
    //////////////////////////////////////////////////////////////////////////
    ////Data Analysis
    ////AnalysisData::writeAnalysisData(para, t);
    //AnalysisData::writeAnalysisDataX(para, t);
    //AnalysisData::writeAnalysisDataZ(para, t);
    //////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //pressure difference
    ////////////////////////////////////////////////////////////////////////
    //if (para->getMyID() == para->getPressInID())       calcPressure(para,  "in", 0);
    //else if (para->getMyID() == para->getPressOutID()) calcPressure(para, "out", 0);
    ////////////////////////////////////////////////////////////////////////
    //flow rate
    ////////////////////////////////////////////////////////////////////////
    //calcFlowRate(para, 0);
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //calculate 2nd, 3rd and higher order moments
    ////////////////////////////////////////////////////////////////////////
    if (para->getCalc2ndOrderMoments())  calc2ndMoments(para.get(), cudaMemoryManager.get());
    if (para->getCalc3rdOrderMoments())  calc3rdMoments(para.get(), cudaMemoryManager.get());
    if (para->getCalcHighOrderMoments()) calcHigherOrderMoments(para.get(), cudaMemoryManager.get());
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //calculate mean on host
    ////////////////////////////////////////////////////////////////////////
    if (para->getCalcMean() && ((int)timestep > para->getTimeCalcMedStart()) && ((int)timestep <= para->getTimeCalcMedEnd()) && ((timestep%(unsigned int)para->getclockCycleForMP())==0))
    {
        unsigned int tdiff = timestep - previousTimestepForAveraging;
        calcMean(para.get(), tdiff);
        /////////////////////////////////
        //added for incremental averaging
        previousTimestepForAveraging = timestep;
        resetMean(para.get());
        /////////////////////////////////
    }
    if (para->getCalcTurbulenceIntensity())
    {
        uint t_diff = timestep - previousTimestepForTurbulenceIntensityCalculation;
        calcTurbulenceIntensity(para.get(), cudaMemoryManager.get(), t_diff);
        //writeAllTiDatafToFile(para.get(), t);
    }
    ////////////////////////////////////////////////////////////////////////
    dataWriter->writeTimestep(para, timestep);
    ////////////////////////////////////////////////////////////////////////
    if (para->getCalcTurbulenceIntensity()) {
        previousTimestepForTurbulenceIntensityCalculation = timestep;
        resetVelocityFluctuationsAndMeans(para.get(), cudaMemoryManager.get());
    }
    ////////////////////////////////////////////////////////////////////////
    if (para->getCalcDragLift()) 
    {
        printDragLift(para.get(), cudaMemoryManager.get(), timestep);
    }
    ////////////////////////////////////////////////////////////////////////

#if DEBUG_FS
    // Write distributions (f's) for debugging purposes.
    DistributionDebugWriter::copyDistributionsToHost(*para, *cudaMemoryManager);
    DistributionDebugWriter::writeDistributions(*para, timestep);
#endif

    ////////////////////////////////////////////////////////////////////////
    VF_LOG_INFO("... done");
    ////////////////////////////////////////////////////////////////////////
}

Simulation::~Simulation()
{
    // Cuda Streams
    if (para->getUseStreams()) {
        para->getStreamManager()->destroyCudaEvents();
        para->getStreamManager()->terminateStreams();
    }

    //CudaFreeHostMemory
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //para->cudaFreeFull(lev);
        cudaMemoryManager->cudaFreeCoord(lev);
        cudaMemoryManager->cudaFreeSP(lev);
        if (para->getCalcMean())
        {
            cudaMemoryManager->cudaFreeMeanSP(lev);
        }
        //para->cudaFreeVeloBC(lev);
        //para->cudaFreeWallBC(lev);
        //para->cudaFreeVeloBC(lev);
        //para->cudaFreeInlet(lev);
        //para->cudaFreeOutlet(lev);
        //para->cudaFreeGeomBC(lev);
        //para->cudaFreePress(lev);
    }
    if (para->getMaxLevel()>1)
    {
        for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
        {
            cudaMemoryManager->cudaFreeInterfaceCF(lev);
            cudaMemoryManager->cudaFreeInterfaceFC(lev);
            cudaMemoryManager->cudaFreeInterfaceOffCF(lev);
            cudaMemoryManager->cudaFreeInterfaceOffFC(lev);
            //para->cudaFreePressX1(lev);
        }
    }
    //para->cudaFreeVeloBC(0); //level = 0
    //para->cudaFreePressBC();
    //para->cudaFreeVeloPropeller(para->getFine());
    //para->cudaFreePressX0(para->getCoarse());

    //////////////////////////////////////////////////////////////////////////
    //Temp
    if (para->getDiffOn() == true)
    {
        for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
        {
            checkCudaErrors(cudaFreeHost(para->getParH(lev)->Conc_Full));
            checkCudaErrors(cudaFreeHost(para->getParH(lev)->concentration));
            checkCudaErrors(cudaFreeHost(para->getParH(lev)->AdvectionDiffusionNoSlipBC.concentration));
            checkCudaErrors(cudaFreeHost(para->getParH(lev)->AdvectionDiffusionNoSlipBC.k));
            checkCudaErrors(cudaFreeHost(para->getParH(lev)->AdvectionDiffusionDirichletBC.concentration));
            checkCudaErrors(cudaFreeHost(para->getParH(lev)->AdvectionDiffusionDirichletBC.k));
        }
    }
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    //free second order moments
    if (para->getCalc2ndOrderMoments())
    {
        for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
        {
            cudaMemoryManager->cudaFree2ndMoments(lev);
        }
    }
    //////////////////////////////////////////////////////////////////////////
    //free third order moments
    if (para->getCalc3rdOrderMoments())
    {
        for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
        {
            cudaMemoryManager->cudaFree3rdMoments(lev);
        }
    }
    //////////////////////////////////////////////////////////////////////////
    //free higher order moments
    if (para->getCalcHighOrderMoments())
    {
        for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
        {
            cudaMemoryManager->cudaFreeHigherMoments(lev);
        }
    }
    //////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////
    //Multi GPU
    //////////////////////////////////////////////////////////////////////////
    ////1D domain decomposition
    //if (para->getNumprocs() > 1)
    //{
    // for (int lev=para->getCoarse(); lev < para->getFine(); lev++)
    // {
    //  for (unsigned int i=0; i < para->getNumberOfProcessNeighbors(lev, "send"); i++)
    //  {
    //   para->cudaFreeProcessNeighbor(lev, i);
    //  }
    // }
    //}
    //////////////////////////////////////////////////////////////////////////
    //3D domain decomposition
    if (para->getNumprocs() > 1)
    {
        for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
        {
            //////////////////////////////////////////////////////////////////////////
            for (unsigned int i = 0; i < para->getNumberOfProcessNeighborsX(lev, "send"); i++)
            {
                cudaMemoryManager->cudaFreeProcessNeighborX(lev, i);
            }
            //////////////////////////////////////////////////////////////////////////
            for (unsigned int i = 0; i < para->getNumberOfProcessNeighborsY(lev, "send"); i++)
            {
                cudaMemoryManager->cudaFreeProcessNeighborY(lev, i);
            }
            //////////////////////////////////////////////////////////////////////////
            for (unsigned int i = 0; i < para->getNumberOfProcessNeighborsZ(lev, "send"); i++)
            {
                cudaMemoryManager->cudaFreeProcessNeighborZ(lev, i);
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////
    //Normals
    if (para->getIsGeoNormal()) {
        for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
        {
            cudaMemoryManager->cudaFreeGeomNormals(lev);
        }
    }
    //////////////////////////////////////////////////////////////////////////
    // Turbulence Intensity
    if (para->getCalcTurbulenceIntensity()) {
        cudaFreeTurbulenceIntensityArrays(para.get(), cudaMemoryManager.get());
    //PreCollisionInteractors
    for( SPtr<PreCollisionInteractor> actuator: para->getActuators()){
        actuator->free(para.get(), cudaMemoryManager.get());
    }

    for( SPtr<PreCollisionInteractor> probe: para->getProbes()){
        probe->free(para.get(), cudaMemoryManager.get());
    }
    //////////////////////////////////////////////////////////////////////////
    }
}
