#include "PhiAndNuTest.h"

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities\TestSimulation\TestSimulation.h"
#include "Utilities\SimulationInfo\SimulationInfo.h"

#include "Tests\PhiAndNuTest\PostProcessingStrategy\PostProcessingStrategyPhiAndNuTest.h"
#include "Tests\PhiAndNuTest\PhiAndNuTestParameterStruct.h"

#include <iomanip>

std::shared_ptr<PhiAndNuTest> PhiAndNuTest::getNewInstance(std::shared_ptr< ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNuTestParameterStruct> testPara, std::string dataToCalculate)
{
	return std::shared_ptr<PhiAndNuTest>(new PhiAndNuTest(colorOutput, viscosity, testPara, dataToCalculate));
}

void PhiAndNuTest::evaluate()
{
	for (int i = 0; i < phiAndNuPostProStrategies.size(); i++) {
		lx.push_back(phiAndNuPostProStrategies.at(i)->getNumberOfXNodes());
		if (dataToCalculate == "Vx") {
			nu.push_back(phiAndNuPostProStrategies.at(i)->getNuVx());
			phiDiff.push_back(phiAndNuPostProStrategies.at(i)->getPhiDiffVx());
		}
		if (dataToCalculate == "Vy") {
			nu.push_back(phiAndNuPostProStrategies.at(i)->getNuVy());
			phiDiff.push_back(phiAndNuPostProStrategies.at(i)->getPhiDiffVy());
		}
		if (dataToCalculate == "Vz") {
			nu.push_back(phiAndNuPostProStrategies.at(i)->getNuVz());
			phiDiff.push_back(phiAndNuPostProStrategies.at(i)->getPhiDiffVz());
		}
	}
	nuDiff = calcNuDiff(nu);
	orderOfAccuracyPhiDiff = calcOrderOfAccuracy(phiDiff);
	orderOfAccuracyNuDiff = calcOrderOfAccuracy(nuDiff);
	phiDiffTestPassed = checkTestPassed(orderOfAccuracyPhiDiff);
	nuDiffTestPassed = checkTestPassed(orderOfAccuracyNuDiff);

	makeConsoleOutput();
}

void PhiAndNuTest::update()
{
	TestImp::update();
}

void PhiAndNuTest::addSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< PhiAndNuTestPostProcessingStrategy> postProStrategy)
{
	TestImp::addSimulation(sim, simInfo, postProStrategy);
	phiAndNuPostProStrategies.push_back(postProStrategy);
}

std::vector<bool> PhiAndNuTest::getPassedTests()
{
	std::vector< bool> passed;
	passed.push_back(phiDiffTestPassed);
	passed.push_back(nuDiffTestPassed);
	return passed;
}

void PhiAndNuTest::makeConsoleOutput()
{
	colorOutput->makeNuTestOutput(nuDiffTestPassed, simInfos.at(0), simInfos.at(1), startStepCalculation, endStepCalculation, dataToCalculate, nu.at(0), nu.at(1), nuDiff.at(0), nuDiff.at(1), orderOfAccuracyNuDiff);
	colorOutput->makePhiTestOutput(phiDiffTestPassed, simInfos.at(0), simInfos.at(1), startStepCalculation, endStepCalculation, dataToCalculate, phiDiff.at(0), phiDiff.at(1), orderOfAccuracyPhiDiff);
}

std::string PhiAndNuTest::getDataToCalculate()
{
	return dataToCalculate;
}

std::vector<int> PhiAndNuTest::getLx()
{
	std::vector<int> lxINT;
	for (int i = 0; i < lx.size(); i++)
		lxINT.push_back((int)lx.at(i));
	return lxINT;
}

std::vector<double> PhiAndNuTest::getNu()
{
	return nu;
}

std::vector<double> PhiAndNuTest::getNuDiff()
{
	return nuDiff;
}

std::vector<double> PhiAndNuTest::getPhiDiff()
{
	return phiDiff;
}

double PhiAndNuTest::getOrderOfAccuracyNuDiff()
{
	return orderOfAccuracyNuDiff;
}

double PhiAndNuTest::getOrderOfAccuracyPhiDiff()
{
	return orderOfAccuracyPhiDiff;
}

PhiAndNuTest::PhiAndNuTest(std::shared_ptr< ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiAndNuTestParameterStruct> testPara, std::string dataToCalculate)
	: TestImp(colorOutput), viscosity(viscosity), dataToCalculate(dataToCalculate)
{
	minOrderOfAccuracy = testPara->minOrderOfAccuracy;
	startStepCalculation = testPara->startTimeStepCalculation;
	endStepCalculation = testPara->endTimeStepCalculation;

	lx.resize(0);
	phiDiff.resize(0);
	nuDiff.resize(0);
}

double PhiAndNuTest::calcOrderOfAccuracy(std::vector<double> data)
{
	double ooa = log(data.at(0) / data.at(1)) / log(lx.at(1) / lx.at(0));
	
	return ooa;
}

bool PhiAndNuTest::checkTestPassed(double orderOfAccuracy)
{
	return orderOfAccuracy > minOrderOfAccuracy;
}

std::vector<double> PhiAndNuTest::calcNuDiff(std::vector<double> nu)
{
	std::vector< double> results;
	for(int i = 0; i < nu.size(); i++)
		results.push_back((nu.at(i) - viscosity) / viscosity);
	return results;
}
