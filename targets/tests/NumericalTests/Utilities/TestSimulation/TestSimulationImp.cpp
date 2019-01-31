#include "TestSimulationImp.h"

#include "Utilities\DataWriter\ToVectorWriter.h"
#include "Utilities\SimulationParameter\SimulationParameter.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Test\SimulationObserver.h"
#include "Utilities\ColorConsoleOutput\ColorConsoleOutput.h"
#include "Utilities\KernelConfiguration\KernelConfiguration.h"
#include "Utilities\DataWriter\AnalyticalResults2DToVTKWriter\AnalyticalResults2DToVTKWriter.h"
#include "Utilities\Structs\TestSimulationDataStruct.h"


#include <sstream>
#include <iomanip>

std::shared_ptr<TestSimulationImp> TestSimulationImp::getNewInsance(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResults, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput)
{
	return std::shared_ptr<TestSimulationImp>(new TestSimulationImp(testSimData, simResults, toVectorWriter, anaResultWriter, colorOutput));
}

TestSimulationImp::TestSimulationImp(std::shared_ptr<TestSimulationDataStruct> testSimData, std::shared_ptr<SimulationResults> simResults, std::shared_ptr<ToVectorWriter> toVectorWriter, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, std::shared_ptr<ColorConsoleOutput> colorOutput)
{
	this->simPara = testSimData->simParameter;
	this->simInfo = testSimData->simInformation;
	this->analyticalResult = testSimData->analyticalResult;
	this->simResults = simResults;
	this->anaResultWriter = anaResultWriter;
	this->toVectorWriter = toVectorWriter;
	this->colorOutput = colorOutput;

	this->simObserver.resize(0);
	this->simualtionRun = false;
}

std::shared_ptr<SimulationParameter> TestSimulationImp::getSimulationParameter()
{
	return simPara;
}

std::shared_ptr<SimulationResults> TestSimulationImp::getSimulationResults()
{
	return simResults;
}

std::shared_ptr<AnalyticalResults> TestSimulationImp::getAnalyticalResults()
{
	return analyticalResult;
}

std::shared_ptr<DataWriter> TestSimulationImp::getDataWriter()
{
	return toVectorWriter;
}

std::shared_ptr<SimulationInfo> TestSimulationImp::getSimulationInfo()
{
	return simInfo;
}

bool TestSimulationImp::getSimulationRun()
{
	return simualtionRun;
}

void TestSimulationImp::registerSimulationObserver(std::shared_ptr< SimulationObserver> simObserver)
{
	this->simObserver.push_back(simObserver);
}

void TestSimulationImp::notifyObserver()
{
	for (int i = 0; i < simObserver.size(); i++)
		simObserver.at(i)->update();
}

float TestSimulationImp::calcTestTime()
{
	float timeInMiliSec = ((float)(testEndTime - testStartTime) / CLOCKS_PER_SEC);
	return timeInMiliSec;
}

void TestSimulationImp::writeAnalyticalResultsToVTK()
{
	if (!analyticalResult->isCalculated())
		analyticalResult->calc(simResults);

	anaResultWriter->writeAnalyticalResult(para, analyticalResult);
}

void TestSimulationImp::makeSimulationHeadOutput()
{
	colorOutput->makeSimulationHeadOutput(simInfo);
}

void TestSimulationImp::setSimulationStartTime()
{
	simulationStartTime = time(NULL);
}

void TestSimulationImp::setSimulationEndTimeAndNotifyObserver()
{
	simulationEndTime = time(NULL);
	simualtionRun = true;

	setTestStartTime();
	notifyObserver();
	setTestEndTime();

	setAnalyticalResultWriteTimeStartTime();
	writeAnalyticalResultsToVTK();
	setAnalyticalResultWriteEndTime();
}

void TestSimulationImp::setTestStartTime()
{
	testStartTime = clock();
}

void TestSimulationImp::setTestEndTime()
{
	testEndTime = clock();
}

void TestSimulationImp::setAnalyticalResultWriteTimeStartTime()
{
	analyticalResultWriteStartTime = time(NULL);
}

void TestSimulationImp::setAnalyticalResultWriteEndTime()
{
	analyticalResultWriteEndTime = time(NULL);
}

double TestSimulationImp::calcSimTime()
{
	return difftime(simulationEndTime, simulationStartTime);
}

double TestSimulationImp::calcAnalyticalResultWriteTime()
{
	return difftime(analyticalResultWriteEndTime, analyticalResultWriteStartTime);
}

std::string TestSimulationImp::getRunTimeOutput()
{
	std::ostringstream oss;
	oss << std::left << std::setfill(' ') << std::setw(11) << "Simulation" << std::setw(17) << simInfo->getSimulationName() << "\t" << std::right << std::setw(3) << simPara->getLx() << "\t\t" << std::setw(9) << calcSimTime() << " sec" << std::endl;
	oss << std::left << std::setfill(' ') << std::setw(11) << "Test" << std::setw(17) << simInfo->getSimulationName() << "\t" << std::right << std::setw(3) << simPara->getLx() << "\t\t" << std::setw(9) << calcTestTime() << " sec" << std::endl;
	oss << std::left << std::setfill(' ') << std::setw(11) << "FileWriting" << std::setw(17) << simInfo->getSimulationName() << "\t" << std::right << std::setw(3) << simPara->getLx() << "\t\t" << std::setw(9) << calcAnalyticalResultWriteTime() << " sec" << std::endl;
	return oss.str();
}

void TestSimulationImp::setParameter(std::shared_ptr<Parameter> para)
{
	this->para = para;
}