#include "TestConditionImp.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU\Output\FileWriter.h"

#include "Utilities\GridReaderforTesting\GridReaderforTesting.h"
#include "Utilities/SimulationResults/SimulationResults.h"
#include "Utilities\DataWriter\Y2dSliceToResults\Y2dSliceToResults.h"
#include "Utilities\InitialCondition\InitialCondition.h"
#include "Utilities\Calculator\Calculator.h"

#include <sstream>

std::shared_ptr<Parameter> TestConditionImp::getParameter()
{
	return para;
}

std::shared_ptr<GridProvider> TestConditionImp::getGrid()
{
	return grid;
}

std::shared_ptr<DataWriter> TestConditionImp::getDataWriter()
{
	return writeToVector;
}

std::shared_ptr<Calculator> TestConditionImp::getCalculator()
{
	return calculator;
}

std::shared_ptr<TestConditionImp> TestConditionImp::getNewInstance()
{
	return std::shared_ptr<TestConditionImp>(new TestConditionImp());
}

void TestConditionImp::initParameter(real viscosity, std::string aGridPath, std::string filePath, int numberOfGridLevels, unsigned int endTime, unsigned int timeStepLength, std::vector<int> devices, real velocity)
{
	para = Parameter::make();

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
}

void TestConditionImp::initInitialConditions(std::shared_ptr<InitialCondition> initialCondition)
{
	this->initialCondition = initialCondition;
	this->initialCondition->setParameter(para);
}

void TestConditionImp::initGridProvider()
{
	grid = std::shared_ptr<GridProvider>(new GridReaderforTesting(para, initialCondition));
}

void TestConditionImp::initCalculator(std::shared_ptr<Calculator> calc)
{
	this->calculator = calc;
}


void TestConditionImp::setTestResults(std::shared_ptr<TestResults> testResults)
{
	this->testResults = testResults;
}


void TestConditionImp::initDataWriter(unsigned int ySliceForCalculation, unsigned int startTimeCalculation, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, unsigned int startTimeDataWriter)
{
	fileWriter = std::shared_ptr<FileWriter>(new FileWriter());
	writeToVector = std::shared_ptr<ToVectorWriter>(new Y2dSliceToResults(simResults, ySliceForCalculation, startTimeCalculation, endTime, timeStepLength, writeFiles, fileWriter, startTimeDataWriter));
}

void TestConditionImp::initSimulationResults(unsigned int lx, unsigned int lz, unsigned int timeStepLength)
{
	simResults = SimulationResults::getNewInstance(lx, lz, timeStepLength);
	calculator->setSimulationResults(simResults);
}
