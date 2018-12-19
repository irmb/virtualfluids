#include "TestSimulationImp.h"

#include "VirtualFluids_GPU\Output\FileWriter.h"

#include "Utilities\SimulationParameter\SimulationParameter.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities\Results\AnalyticalResults\AnalyticalResult.h"
#include "Utilities\Test\SimulationObserver.h"
#include "Utilities\DataWriter\Y2dSliceToResults\Y2dSliceToResults.h"
#include "Utilities\ColorConsoleOutput\ColorConsoleOutput.h"
#include "Utilities\KernelConfiguration\KernelConfiguration.h"
#include "Utilities\DataWriter\AnalyticalResults2DToVTKWriter\AnalyticalResults2DToVTKWriter.h"

#include <sstream>
#include <iomanip>

std::shared_ptr<TestSimulation> TestSimulationImp::getNewInsance(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< ColorConsoleOutput> colorOutput, std::shared_ptr< SimulationResults> simResults, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, bool writeAnalyticalResults)
{
	return std::shared_ptr< TestSimulation>(new TestSimulationImp(simID, simPara, simInfo, colorOutput, simResults, analyticalResult, anaResultWriter, writeAnalyticalResults));
}

std::shared_ptr<SimulationParameter> TestSimulationImp::getSimulationParameter()
{
	return simPara;
}

std::shared_ptr<SimulationResults> TestSimulationImp::getSimulationResults()
{
	return simResults;
}

std::shared_ptr<DataWriter> TestSimulationImp::getDataWriter()
{
	return writeToVector;
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

	if (writeAnalyticalResults) {
		setAnalyticalResultWriteTimeStartTime();
		writeAnalyticalResultsToVTK();
		setAnalyticalResultWriteEndTime();
	}
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
	if (writeAnalyticalResults)
		oss << std::left << std::setfill(' ') << std::setw(11) << "FileWriting" << std::setw(17) << simInfo->getSimulationName() << "\t" << std::right << std::setw(3) << simPara->getLx() << "\t\t" << std::setw(9) << calcAnalyticalResultWriteTime() << " sec" << std::endl;
	return oss.str();
}

void TestSimulationImp::setParameter(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

TestSimulationImp::TestSimulationImp(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< ColorConsoleOutput> colorOutput, std::shared_ptr< SimulationResults> simResults, std::shared_ptr<AnalyticalResults> analyticalResult, std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter, bool writeAnalyticalResults) : simID(simID), colorOutput(colorOutput), simResults(simResults), writeAnalyticalResults(writeAnalyticalResults)
{
	this->simPara = simPara;
	this->simInfo = simInfo;
	this->simInfo->setSimulationID(simID);
	this->analyticalResult = analyticalResult;
	this->anaResultWriter = anaResultWriter;

	
	writeToVector = Y2dSliceToResults::getNewInstance(simResults, simPara->getYSliceForCalculation(), simPara->getStartTimeCalculation(), simPara->getEndTime(), simPara->getTimeStepLength(), simPara->getWriteFiles(), std::shared_ptr<FileWriter>(new FileWriter()), simPara->getStartTimeDataWriter());

	simObserver.resize(0);
	simualtionRun = false;
}