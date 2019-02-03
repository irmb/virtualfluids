#include "PhiAndNyTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\TestSimulation\TestSimulation.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"

#include "Tests\PhiAndNyTest\PostProcessingStrategy\PostProcessingStrategyPhiAndNyTest.h"
#include "Tests\PhiAndNyTest\PhiAndNyTestParameterStruct.h"

#include <iomanip>

std::shared_ptr<PhiAndNyTest> PhiAndNyTest::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNyTestParameterStruct> testPara, std::string dataToCalculate)
{
	return std::shared_ptr<PhiAndNyTest>(new PhiAndNyTest(colorOutput, viscosity, testPara, dataToCalculate));
}

void PhiAndNyTest::evaluate()
{
	for (int i = 0; i < phiAndNyPostProStrategies.size(); i++) {
		lx.push_back(phiAndNyPostProStrategies.at(i)->getNumberOfXNodes());
		if (dataToCalculate == "Vx") {
			ny.push_back(phiAndNyPostProStrategies.at(i)->getNyVx());
			phiDiff.push_back(phiAndNyPostProStrategies.at(i)->getPhiDiffVx());
		}
		if (dataToCalculate == "Vy") {
			ny.push_back(phiAndNyPostProStrategies.at(i)->getNyVy());
			phiDiff.push_back(phiAndNyPostProStrategies.at(i)->getPhiDiffVy());
		}
		if (dataToCalculate == "Vz") {
			ny.push_back(phiAndNyPostProStrategies.at(i)->getNyVz());
			phiDiff.push_back(phiAndNyPostProStrategies.at(i)->getPhiDiffVz());
		}
	}
	nyDiff = calcNyDiff(ny);
	orderOfAccuracyPhiDiff = calcOrderOfAccuracy(phiDiff);
	orderOfAccuracyNyDiff = calcOrderOfAccuracy(nyDiff);
	phiDiffTestPassed = checkTestPassed(orderOfAccuracyPhiDiff);
	nyDiffTestPassed = checkTestPassed(orderOfAccuracyNyDiff);

	makeConsoleOutput();
}

void PhiAndNyTest::update()
{
	TestImp::update();
}

void PhiAndNyTest::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PhiAndNyTestPostProcessingStrategy> postProStrategy)
{
	TestImp::addSimulation(sim, simInfo, postProStrategy);
	phiAndNyPostProStrategies.push_back(postProStrategy);
}

std::vector<bool> PhiAndNyTest::getPassedTests()
{
	std::vector<bool> passed;
	passed.push_back(phiDiffTestPassed);
	passed.push_back(nyDiffTestPassed);
	return passed;
}

void PhiAndNyTest::makeConsoleOutput()
{
	colorOutput->makeNyTestOutput(nyDiffTestPassed, simInfos.at(0), simInfos.at(1), startStepCalculation, endStepCalculation, dataToCalculate, ny.at(0), ny.at(1), nyDiff.at(0), nyDiff.at(1), orderOfAccuracyNyDiff);
	colorOutput->makePhiTestOutput(phiDiffTestPassed, simInfos.at(0), simInfos.at(1), startStepCalculation, endStepCalculation, dataToCalculate, phiDiff.at(0), phiDiff.at(1), orderOfAccuracyPhiDiff);
}

std::string PhiAndNyTest::getDataToCalculate()
{
	return dataToCalculate;
}

std::vector<int> PhiAndNyTest::getLx()
{
	std::vector<int> lxINT;
	for (int i = 0; i < lx.size(); i++)
		lxINT.push_back((int)lx.at(i));
	return lxINT;
}

std::vector<double> PhiAndNyTest::getNy()
{
	return ny;
}

std::vector<double> PhiAndNyTest::getNyDiff()
{
	return nyDiff;
}

std::vector<double> PhiAndNyTest::getPhiDiff()
{
	return phiDiff;
}

double PhiAndNyTest::getOrderOfAccuracyNyDiff()
{
	return orderOfAccuracyNyDiff;
}

double PhiAndNyTest::getOrderOfAccuracyPhiDiff()
{
	return orderOfAccuracyPhiDiff;
}

PhiAndNyTest::PhiAndNyTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNyTestParameterStruct> testPara, std::string dataToCalculate)
	: TestImp(colorOutput), viscosity(viscosity), dataToCalculate(dataToCalculate)
{
	minOrderOfAccuracy = testPara->minOrderOfAccuracy;
	startStepCalculation = testPara->startTimeStepCalculation;
	endStepCalculation = testPara->endTimeStepCalculation;

	lx.resize(0);
	phiDiff.resize(0);
	nyDiff.resize(0);
}

double PhiAndNyTest::calcOrderOfAccuracy(std::vector<double> data)
{
	double ooa = log(data.at(0) / data.at(1)) / log(lx.at(1) / lx.at(0));
	
	return ooa;
}

bool PhiAndNyTest::checkTestPassed(double orderOfAccuracy)
{
	return orderOfAccuracy > minOrderOfAccuracy;
}

std::vector<double> PhiAndNyTest::calcNyDiff(std::vector<double> ny)
{
	std::vector<double> results;
	for(int i = 0; i < ny.size(); i++)
		results.push_back((ny.at(i) - viscosity) / viscosity);
	return results;
}
