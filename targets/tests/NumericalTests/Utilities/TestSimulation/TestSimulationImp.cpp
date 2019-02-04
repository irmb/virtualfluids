#include "TestSimulationImp.h"

#include "Utilities\DataWriter\ToVectorWriter.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Test\SimulationObserver.h"
#include "Utilities\ColorConsoleOutput\ColorConsoleOutput.h"
#include "Utilities\KernelConfiguration\KernelConfiguration.h"
#include "Utilities\DataWriter\AnalyticalResults2DToVTKWriter\AnalyticalResults2DToVTKWriter.h"
#include "Utilities\Structs\TestSimulationDataStruct.h"
#include "Utilities\Time\TimeTracking.h"


std::shared_ptr<TestSimulationImp> TestSimulationImp::getNewInsance(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResult, std::shared_ptr<TimeTracking> timeTracking, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput)
{
	return std::shared_ptr<TestSimulationImp>(new TestSimulationImp(testSimData, simResult, timeTracking, toVectorWriter, anaResultWriter, colorOutput));
}

TestSimulationImp::TestSimulationImp(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResult, std::shared_ptr<TimeTracking> timeTracking, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput)
{
	this->simPara = testSimData->simParameter;
	this->simInfo = testSimData->simInformation;
	this->analyticalResult = testSimData->analyticalResult;
	this->initialCondition = testSimData->initialCondition;

	this->timeTracking = timeTracking;

	this->simResult = simResult;
	this->toVectorWriter = toVectorWriter;

	this->anaResultWriter = anaResultWriter;
	this->colorOutput = colorOutput;
	
	this->simObserver.resize(0);
	this->simualtionRun = false;
}

std::shared_ptr<SimulationParameter> TestSimulationImp::getSimulationParameter()
{
	return simPara;
}

std::shared_ptr<AnalyticalResults> TestSimulationImp::getAnalyticalResults()
{
	return analyticalResult;
}

std::shared_ptr<DataWriter> TestSimulationImp::getDataWriter()
{
	return toVectorWriter;
}

std::shared_ptr<InitialCondition> TestSimulationImp::getInitialCondition()
{
	return initialCondition;
}

std::shared_ptr<SimulationInfo> TestSimulationImp::getSimulationInfo()
{
	return simInfo;
}

std::shared_ptr<TimeTracking> TestSimulationImp::getTimeTracking()
{
	return timeTracking;
}

bool TestSimulationImp::getSimulationRun()
{
	return simualtionRun;
}

void TestSimulationImp::registerSimulationObserver(std::shared_ptr<SimulationObserver> simObserver)
{
	this->simObserver.push_back(simObserver);
}

void TestSimulationImp::notifyObserver()
{
	for (int i = 0; i < simObserver.size(); i++)
		simObserver.at(i)->update();
}

void TestSimulationImp::writeAnalyticalResultsToVTK()
{
	if (!analyticalResult->isCalculated())
		analyticalResult->calc(simResult);

	anaResultWriter->writeAnalyticalResult(para, analyticalResult);
}

void TestSimulationImp::makeSimulationHeadOutput()
{
	colorOutput->makeSimulationHeadOutput(simInfo);
}

void TestSimulationImp::startPostProcessing()
{
	simualtionRun = true;

	timeTracking->setTestStartTime();
	notifyObserver();
	timeTracking->setTestEndTime();

	timeTracking->setAnalyticalResultWriteStartTime();
	writeAnalyticalResultsToVTK();
	timeTracking->setAnalyticalResultWriteEndTime();
}

void TestSimulationImp::setParameter(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

std::shared_ptr<SimulationResults> TestSimulationImp::getSimulationResults()
{
	return simResult;
}