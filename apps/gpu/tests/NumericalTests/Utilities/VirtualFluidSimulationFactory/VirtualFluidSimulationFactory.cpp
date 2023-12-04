#include "VirtualFluidSimulationFactory.h"

#include "Utilities/InitialCondition/InitialCondition.h"
#include "Utilities/KernelConfiguration/KernelConfiguration.h"
#include "Utilities/NumericalTestGridReader/NumericalTestGridReader.h"
#include "Utilities/SimulationParameter/SimulationParameter.h"

#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/Parameter/Parameter.h"

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/LBM/Simulation.h"

#include <parallel/MPICommunicator.h>

std::shared_ptr<Parameter> vf::gpu::tests::makeParameter(std::shared_ptr<SimulationParameter> simPara)
{
    auto para = std::make_shared<Parameter>(1, 0);

    para->setQuadricLimiters(0.01, 0.01, 0.01);

    para->setMaxDev(simPara->getDevices().size());
    para->setDevices(simPara->getDevices());

    std::string _prefix = "cells";
    std::string gridPath = simPara->getGridPath() + "/";
    para->setOutputPath(simPara->getFilePath());
    para->setOutputPrefix(_prefix);
    para->setPrintFiles(true);

    para->setD3Qxx(27);
    para->setMaxLevel(simPara->getNumberOfGridLevels());

    para->setTimestepEnd(simPara->getEndTime());
    para->setTimestepOut(simPara->getTimeStepLength());
    para->setTimestepStartOut(1);

    para->setViscosityLB(simPara->getViscosity());
    para->setVelocityLB(simPara->getMaxVelocity());
    para->setViscosityRatio(1.0);
    para->setVelocityRatio(1.0);
    para->setDensityRatio(1.0);
    para->setFactorPressBC(100000.0);

    para->setgeoVec(gridPath + "geoVec.dat");
    para->setcoordX(gridPath + "coordX.dat");
    para->setcoordY(gridPath + "coordY.dat");
    para->setcoordZ(gridPath + "coordZ.dat");
    para->setneighborX(gridPath + "neighborX.dat");
    para->setneighborY(gridPath + "neighborY.dat");
    para->setneighborZ(gridPath + "neighborZ.dat");
    para->setneighborWSB(gridPath + "neighborWSB.dat");
    para->setgeomBoundaryBcQs(gridPath + "geomBoundaryQs.dat");
    para->setgeomBoundaryBcValues(gridPath + "geomBoundaryValues.dat");
    para->setinletBcQs(gridPath + "inletBoundaryQs.dat");
    para->setinletBcValues(gridPath + "inletBoundaryValues.dat");
    para->setoutletBcQs(gridPath + "outletBoundaryQs.dat");
    para->setoutletBcValues(gridPath + "outletBoundaryValues.dat");
    para->settopBcQs(gridPath + "topBoundaryQs.dat");
    para->settopBcValues(gridPath + "topBoundaryValues.dat");
    para->setbottomBcQs(gridPath + "bottomBoundaryQs.dat");
    para->setbottomBcValues(gridPath + "bottomBoundaryValues.dat");
    para->setfrontBcQs(gridPath + "frontBoundaryQs.dat");
    para->setfrontBcValues(gridPath + "frontBoundaryValues.dat");
    para->setbackBcQs(gridPath + "backBoundaryQs.dat");
    para->setbackBcValues(gridPath + "backBoundaryValues.dat");
    para->setnumberNodes(gridPath + "numberNodes.dat");
    para->setLBMvsSI(gridPath + "LBMvsSI.dat");
    para->setscaleCFC(gridPath + "scaleCFC.dat");
    para->setscaleCFF(gridPath + "scaleCFF.dat");
    para->setscaleFCC(gridPath + "scaleFCC.dat");
    para->setscaleFCF(gridPath + "scaleFCF.dat");
    para->setscaleOffsetCF(gridPath + "offsetVecCF.dat");
    para->setscaleOffsetFC(gridPath + "offsetVecFC.dat");
    para->setDiffOn(false);
    para->setDoCheckPoint(false);
    para->setDoRestart(false);
    para->setUseGeometryValues(false);
    para->setCalc2ndOrderMoments(false);
    para->setCalc3rdOrderMoments(false);
    para->setCalcHighOrderMoments(false);
    para->setReadGeo(false);
    para->setCalcMean(false);
    para->setConcFile(false);
    para->setUseMeasurePoints(false);
    para->setForcing(0.0, 0.0, 0.0);

    std::vector<int> dist;
    dist.resize(1);
    dist[0] = 0;
    para->setDistX(dist);
    para->setDistY(dist);
    para->setDistZ(dist);

    // TODO: Find out if we still need this
    // para->setNeedInterface(std::vector<bool>{true, true, true, true, true, true});

    para->configureMainKernel(simPara->getKernelConfiguration()->getMainKernel());
    para->setMultiKernelOn(simPara->getKernelConfiguration()->getMultiKernelOn());
    para->setMultiKernelLevel(simPara->getKernelConfiguration()->getMultiKernelLevel());
    para->setMultiKernel(simPara->getKernelConfiguration()->getMultiKernel());

    return para;
}

std::shared_ptr<NumericalTestGridReader> makeGridReader(std::shared_ptr<InitialCondition> initialCondition,
                                                        std::shared_ptr<Parameter> para,
                                                        std::shared_ptr<CudaMemoryManager> cudaManager)
{
    return NumericalTestGridReader::getNewInstance(para, initialCondition, cudaManager);
}

const std::function<void()> vf::gpu::tests::makeVirtualFluidSimulation(std::shared_ptr<Parameter> para,
                                                                       std::shared_ptr<InitialCondition> condition,
                                                                       std::shared_ptr<DataWriter> dataWriter)
{
    auto cudaManager = std::make_shared<CudaMemoryManager>(para);
    auto grid = makeGridReader(condition, para, cudaManager);
    BoundaryConditionFactory bc_factory;
    vf::parallel::Communicator &communicator = *vf::parallel::MPICommunicator::getInstance();
    auto simulation =
        std::make_shared<Simulation>(para, cudaManager, communicator, *grid.get(), &bc_factory);
    simulation->setDataWriter(dataWriter);

    return [simulation]() { simulation->run(); };
}
