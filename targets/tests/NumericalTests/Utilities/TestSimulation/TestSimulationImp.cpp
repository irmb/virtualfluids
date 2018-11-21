#include "TestSimulationImp.h"

#include "VirtualFluids_GPU\Output\FileWriter.h"

#include "Utilities\SimulationParameter\SimulationParameter.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"
#include "Utilities/Results/SimulationResults/SimulationResults.h"
#include "Utilities\Test\SimulationObserver.h"
#include "Utilities\DataWriter\Y2dSliceToResults\Y2dSliceToResults.h"
#include "Utilities\TestCout\TestCoutImp.h"
#include "Utilities\KernelConfiguration\KernelConfiguration.h"

#include <sstream>
#include <iomanip>


std::shared_ptr<TestSimulation> TestSimulationImp::getNewInsance(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo)
{
	return std::shared_ptr< TestSimulation>(new TestSimulationImp(simID, simPara, simInfo));
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

double TestSimulationImp::calcSimTime()
{
	return difftime(endTime, startTime);
}

void TestSimulationImp::makeSimulationHeadOutput()
{
	colorOutput->makeSimulationHeadOutput(simInfo);
}

void TestSimulationImp::setStartTime()
{
	startTime = time(NULL);
}

void TestSimulationImp::setEndTime()
{
	endTime = time(NULL);
	simualtionRun = true;
	notifyObserver();
}

std::string TestSimulationImp::getSimulationRunTimeOutput()
{
	std::ostringstream oss;
	oss << std::left << std::setfill(' ') << std::setw(17) << simInfo->getSimulationName() << "\t" << std::right << std::setw(3) << simPara->getLx() << "\t\t" << std::setw(9) << calcSimTime() << " sec" << std::endl;
	return oss.str();
}

TestSimulationImp::TestSimulationImp(int simID, std::shared_ptr< SimulationParameter> simPara, std::shared_ptr< SimulationInfo> simInfo) : simID(simID)
{
	this->simPara = simPara;
	this->simInfo = simInfo;
	simResults = SimulationResults::getNewInstance(simPara->getLx(), 1, simPara->getLz(), simPara->getTimeStepLength());
	
	writeToVector = std::shared_ptr<ToVectorWriter>(new Y2dSliceToResults(simResults, simPara->getYSliceForCalculation(), simPara->getStartTimeCalculation(), simPara->getEndTime(), simPara->getTimeStepLength(), simPara->getWriteFiles(), std::shared_ptr<FileWriter>(new FileWriter()), simPara->getStartTimeDataWriter()));
	colorOutput = TestCoutImp::getNewInstance();

	simObserver.resize(0);
	simualtionRun = false;
}