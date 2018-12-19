#include "VirtualFluidSimulationImp.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"

#include "Utilities\GridReaderforTesting\GridReaderforTesting.h"
#include "Utilities\InitialCondition\InitialCondition.h"
#include "Utilities\KernelConfiguration\KernelConfiguration.h"

#include "Utilities\TestSimulation\TestSimulation.h"

#include <sstream>

void VirtualFluidSimulationImp::run()
{
	testSim->makeSimulationHeadOutput();
	testSim->setSimulationStartTime();

	Simulation sim;
	sim.init(para, grid, dataWriter);
	sim.run();

	testSim->setSimulationEndTimeAndNotifyObserver();

	sim.free();
}

std::shared_ptr<VirtualFluidSimulationImp> VirtualFluidSimulationImp::getNewInstance()
{
	return std::shared_ptr<VirtualFluidSimulationImp>(new VirtualFluidSimulationImp());
}

void VirtualFluidSimulationImp::initParameter(std::shared_ptr<Parameter> para, std::shared_ptr< KernelConfiguration> kernelConfig, real viscosity, std::string aGridPath, std::string filePath, int numberOfGridLevels, unsigned int endTime, unsigned int timeStepLength, std::vector<int> devices, real velocity)
{
	this->para = para;

	para->setMaxDev(devices.size());
	para->setDevices(devices);
	para->setNumprocs(1);

	std::string _prefix = "cells";
	std::string gridPath = aGridPath + "\\";
	para->setFName(filePath + "/" + _prefix);
	para->setPrintFiles(true);

	para->setD3Qxx(27);
	para->setMaxLevel(numberOfGridLevels);

	para->setTEnd(endTime);
	para->setTOut(timeStepLength);
    para->setTStartOut(1);

	para->setViscosity(viscosity);
	para->setVelocity(velocity);
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
    para->setCalcParticles(false);
    para->setDiffOn(false);
    para->setDoCheckPoint(false);
    para->setDoRestart(false);
    para->setGeometryValues(false);
    para->setCalc2ndOrderMoments(false);
    para->setCalc3rdOrderMoments(false);
    para->setCalcHighOrderMoments(false);
    para->setReadGeo(false);
    para->setCalcMedian(false);
    para->setConcFile(false);
    para->setUseMeasurePoints(false);
    para->setUseWale(false);
    para->setSimulatePorousMedia(false);
	para->setForcing(0.0, 0.0, 0.0);

	std::vector<int> dist;
	dist.resize(1);
	dist[0] = 0;
	para->setDistX(dist);
	para->setDistY(dist);
	para->setDistZ(dist);

    para->setNeedInterface(std::vector<bool>{true, true, true, true, true, true});

	para->setMainKernel(kernelConfig->getMainKernel());
	para->setMultiKernelOn(kernelConfig->getMultiKernelOn());
	para->setMultiKernelLevel(kernelConfig->getMultiKernelLevel());
	para->setMultiKernelName(kernelConfig->getMultiKernelName());
}

void VirtualFluidSimulationImp::initInitialConditions(std::shared_ptr<InitialCondition> initialCondition)
{
	this->initialCondition = initialCondition;
	this->initialCondition->setParameter(para);
}

void VirtualFluidSimulationImp::initGridProvider()
{
	grid = std::shared_ptr<GridProvider>(new GridReaderforTesting(para, initialCondition));
}

void VirtualFluidSimulationImp::setDataWriter(std::shared_ptr<DataWriter> dataWriter)
{
	this->dataWriter = dataWriter;
}

void VirtualFluidSimulationImp::setTestSimulation(std::shared_ptr<TestSimulation> testSim)
{
	this->testSim = testSim;
}
